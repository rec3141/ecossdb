#!/usr/bin/env python3
"""Parse CICES V5.2 xlsx into standardized TSV + JSON hierarchy.

Reads the official CICES spreadsheet and produces:
  - es_ontology.tsv: flat TSV with explicit hierarchy columns
  - es_hierarchy.json: nested JSON tree for visualization

Usage:
    parse_cices.py --input CICES_V5.2.xlsx --output-tsv es_ontology.tsv --output-json es_hierarchy.json
"""

import argparse
import json
import sys
from collections import OrderedDict

import openpyxl


def parse_cices_xlsx(xlsx_path):
    """Parse CICES xlsx into list of class-level records."""
    wb = openpyxl.load_workbook(xlsx_path, read_only=True)
    ws = wb['CICES V5.2']

    records = []
    header_row = None

    for i, row in enumerate(ws.iter_rows(values_only=True)):
        # Find header row by looking for 'Section' column
        if header_row is None:
            vals = [str(v).strip() if v else '' for v in row]
            if 'Section' in vals and 'Division' in vals:
                header_row = i
                continue
            continue

        # Skip empty rows
        if not row or not any(row[:6]):
            continue

        filter_col = str(row[0] or '').strip()
        section = str(row[1] or '').strip()
        division = str(row[2] or '').strip()
        group = str(row[3] or '').strip()
        class_name = str(row[4] or '').strip()
        code = str(row[5] or '').strip()
        class_type = str(row[6] or '').strip()
        guidance = str(row[7] or '').strip()
        simple_desc = str(row[11] or '').strip() if len(row) > 11 else ''

        if not code:
            continue

        records.append({
            'code': code,
            'filter': filter_col,
            'section': section,
            'division': division,
            'group': group,
            'class': class_name,
            'class_type': class_type,
            'simple_descriptor': simple_desc,
            'guidance': guidance,
        })

    wb.close()
    return records


def build_hierarchy(records):
    """Build full hierarchy from class-level records.

    CICES codes are hierarchical: 1.1.1.1 implies section 1, division 1.1,
    group 1.1.1. We extract each level and build the full ontology.
    """
    # Collect unique hierarchy entries
    entries = OrderedDict()  # code -> entry dict

    for rec in records:
        code = rec['code']
        parts = code.split('.')

        # Clean section name (remove parenthetical qualifiers)
        section_name = rec['section'].split('(')[0].strip()

        # Build each level
        if len(parts) >= 1:
            sec_code = parts[0]
            if sec_code not in entries:
                entries[sec_code] = {
                    'es_code': sec_code,
                    'level': 'section',
                    'section': section_name,
                    'division': '-',
                    'group': '-',
                    'class': '-',
                    'name': section_name,
                    'description': section_name,
                    'parent_code': '-',
                }

        if len(parts) >= 2:
            div_code = '.'.join(parts[:2])
            if div_code not in entries:
                entries[div_code] = {
                    'es_code': div_code,
                    'level': 'division',
                    'section': section_name,
                    'division': rec['division'],
                    'group': '-',
                    'class': '-',
                    'name': rec['division'],
                    'description': rec['division'],
                    'parent_code': parts[0],
                }

        if len(parts) >= 3:
            grp_code = '.'.join(parts[:3])
            if grp_code not in entries:
                entries[grp_code] = {
                    'es_code': grp_code,
                    'level': 'group',
                    'section': section_name,
                    'division': rec['division'],
                    'group': rec['group'],
                    'class': '-',
                    'name': rec['group'],
                    'description': rec['group'],
                    'parent_code': '.'.join(parts[:2]),
                }

        if len(parts) >= 4:
            cls_code = '.'.join(parts[:4])
            if cls_code not in entries:
                desc = rec['simple_descriptor'] or rec['class']
                entries[cls_code] = {
                    'es_code': cls_code,
                    'level': 'class',
                    'section': section_name,
                    'division': rec['division'],
                    'group': rec['group'],
                    'class': rec['class'],
                    'name': desc,
                    'description': rec['class'],
                    'parent_code': '.'.join(parts[:3]),
                }

        # 5th digit extensions (class types)
        if len(parts) >= 5:
            ext_code = code
            if ext_code not in entries:
                entries[ext_code] = {
                    'es_code': ext_code,
                    'level': 'subclass',
                    'section': section_name,
                    'division': rec['division'],
                    'group': rec['group'],
                    'class': rec['class'],
                    'name': rec['class_type'] or rec['class'],
                    'description': rec['class_type'] or rec['class'],
                    'parent_code': '.'.join(parts[:4]),
                }

    return list(entries.values())


def write_tsv(entries, output_path):
    """Write flat TSV ontology."""
    cols = ['es_code', 'level', 'section', 'division', 'group', 'class',
            'name', 'description', 'parent_code']
    with open(output_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for entry in entries:
            f.write('\t'.join(str(entry.get(c, '-')) for c in cols) + '\n')


def build_json_tree(entries):
    """Build nested JSON tree for viz treemap."""
    by_code = {e['es_code']: e for e in entries}
    tree = {'name': 'Ecosystem Services', 'code': 'root', 'children': []}

    nodes = {'root': tree}

    # Sort by code to ensure parents are processed first
    for entry in sorted(entries, key=lambda e: e['es_code']):
        code = entry['es_code']
        parent = entry['parent_code']

        node = {
            'name': entry['name'],
            'code': code,
            'level': entry['level'],
            'children': [],
        }
        nodes[code] = node

        parent_node = nodes.get(parent, tree)
        parent_node['children'].append(node)

    # Remove empty children arrays from leaves
    def prune(node):
        if not node['children']:
            del node['children']
        else:
            for child in node['children']:
                prune(child)

    prune(tree)
    return tree


def main():
    parser = argparse.ArgumentParser(description='Parse CICES V5.2 xlsx')
    parser.add_argument('--input', required=True, help='CICES xlsx file')
    parser.add_argument('--output-tsv', required=True, help='Output TSV')
    parser.add_argument('--output-json', required=True, help='Output JSON tree')
    args = parser.parse_args()

    records = parse_cices_xlsx(args.input)
    if not records:
        print('ERROR: No CICES records found in xlsx', file=sys.stderr)
        sys.exit(1)

    entries = build_hierarchy(records)
    write_tsv(entries, args.output_tsv)

    tree = build_json_tree(entries)
    with open(args.output_json, 'w') as f:
        json.dump(tree, f, indent=2)

    print(f'Parsed {len(records)} CICES classes into {len(entries)} hierarchy entries')


if __name__ == '__main__':
    main()
