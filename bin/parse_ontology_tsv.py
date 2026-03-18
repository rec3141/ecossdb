#!/usr/bin/env python3
"""Parse a generic hierarchical ES ontology TSV into standardized format.

Accepts any TSV with columns: es_code, level, section, division, group, class,
name, description, parent_code — and produces the same format plus a JSON tree.
Also handles simpler formats with just: es_code, name, parent_code, level.

Usage:
    parse_ontology_tsv.py --input custom_ontology.tsv --output-tsv es_ontology.tsv --output-json es_hierarchy.json
"""

import argparse
import csv
import json
import sys


REQUIRED_MINIMAL = {'es_code', 'name', 'parent_code'}
FULL_COLUMNS = ['es_code', 'level', 'section', 'division', 'group', 'class',
                'name', 'description', 'parent_code']


def parse_tsv(input_path):
    """Parse input TSV, handling both full and minimal formats."""
    entries = []
    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        cols = set(reader.fieldnames or [])

        if not REQUIRED_MINIMAL.issubset(cols):
            missing = REQUIRED_MINIMAL - cols
            print(f'ERROR: Missing required columns: {missing}', file=sys.stderr)
            sys.exit(1)

        for row in reader:
            entry = {}
            for col in FULL_COLUMNS:
                entry[col] = row.get(col, '-').strip() or '-'
            entries.append(entry)

    return entries


def infer_levels(entries):
    """Infer level from code depth if not provided."""
    level_map = {1: 'section', 2: 'division', 3: 'group', 4: 'class', 5: 'subclass'}
    for entry in entries:
        if entry['level'] == '-':
            depth = len(entry['es_code'].split('.'))
            entry['level'] = level_map.get(depth, f'level_{depth}')
    return entries


def build_json_tree(entries):
    """Build nested JSON tree."""
    tree = {'name': 'Ecosystem Services', 'code': 'root', 'children': []}
    nodes = {'root': tree}

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

    def prune(node):
        if not node.get('children'):
            node.pop('children', None)
        else:
            for child in node['children']:
                prune(child)

    prune(tree)
    return tree


def main():
    parser = argparse.ArgumentParser(description='Parse generic ES ontology TSV')
    parser.add_argument('--input', required=True, help='Input TSV')
    parser.add_argument('--output-tsv', required=True, help='Output standardized TSV')
    parser.add_argument('--output-json', required=True, help='Output JSON tree')
    args = parser.parse_args()

    entries = parse_tsv(args.input)
    entries = infer_levels(entries)

    with open(args.output_tsv, 'w') as f:
        f.write('\t'.join(FULL_COLUMNS) + '\n')
        for entry in entries:
            f.write('\t'.join(entry[c] for c in FULL_COLUMNS) + '\n')

    tree = build_json_tree(entries)
    with open(args.output_json, 'w') as f:
        json.dump(tree, f, indent=2)

    print(f'Parsed {len(entries)} ontology entries')


if __name__ == '__main__':
    main()
