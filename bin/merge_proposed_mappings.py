#!/usr/bin/env python3
"""Merge agent-proposed mappings into the master mapping table.

Handles:
  - Exact duplicates (skip)
  - Same gene+ES but different confidence (keep higher, log conflict)
  - Same gene+ES but different functional_role (flag for human review)

Usage:
    merge_proposed_mappings.py \
        --master db/mappings/es_gene_mapping.tsv \
        --proposed db/mappings/proposed/agent_*.tsv \
        --output db/mappings/es_gene_mapping_merged.tsv \
        --log merge_log.tsv
"""

import argparse
import csv
import glob
import sys
from collections import defaultdict


COLUMNS = ['es_code', 'es_name', 'gene_id', 'gene_id_type', 'confidence',
           'evidence_source', 'functional_role', 'pathway_context', 'notes']


def load_mapping(filepath):
    """Load mapping TSV into list of row dicts."""
    rows = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def key_of(row):
    return (row.get('gene_id', '').strip(), row.get('es_code', '').strip())


def main():
    parser = argparse.ArgumentParser(description='Merge proposed ES mappings')
    parser.add_argument('--master', required=True, help='Master mapping TSV')
    parser.add_argument('--proposed', required=True, nargs='+',
                        help='Proposed mapping TSV(s) or glob patterns')
    parser.add_argument('--output', required=True, help='Output merged mapping TSV')
    parser.add_argument('--log', required=True, help='Merge log TSV')
    args = parser.parse_args()

    # Load master
    master_rows = load_mapping(args.master)
    master_index = {}
    for row in master_rows:
        master_index[key_of(row)] = row

    # Expand glob patterns
    proposed_files = []
    for pattern in args.proposed:
        expanded = glob.glob(pattern)
        proposed_files.extend(expanded if expanded else [pattern])

    log_entries = []
    added = 0
    skipped = 0
    conflicts = 0

    for filepath in proposed_files:
        proposed_rows = load_mapping(filepath)
        for row in proposed_rows:
            key = key_of(row)

            if key in master_index:
                existing = master_index[key]
                existing_conf = float(existing.get('confidence', 0))
                new_conf = float(row.get('confidence', 0))

                # Check for role conflict
                if existing.get('functional_role') != row.get('functional_role'):
                    conflicts += 1
                    log_entries.append({
                        'action': 'ROLE_CONFLICT',
                        'gene_id': row.get('gene_id', ''),
                        'es_code': row.get('es_code', ''),
                        'existing_role': existing.get('functional_role', ''),
                        'proposed_role': row.get('functional_role', ''),
                        'source': filepath,
                        'note': 'Requires human review',
                    })
                    continue

                # Same role — keep higher confidence
                if new_conf > existing_conf:
                    master_index[key] = row
                    log_entries.append({
                        'action': 'UPDATED_CONFIDENCE',
                        'gene_id': row.get('gene_id', ''),
                        'es_code': row.get('es_code', ''),
                        'old_confidence': str(existing_conf),
                        'new_confidence': str(new_conf),
                        'source': filepath,
                        'note': '',
                    })
                else:
                    skipped += 1
                    log_entries.append({
                        'action': 'SKIPPED_DUPLICATE',
                        'gene_id': row.get('gene_id', ''),
                        'es_code': row.get('es_code', ''),
                        'source': filepath,
                        'note': f'Existing confidence {existing_conf} >= proposed {new_conf}',
                    })
            else:
                master_index[key] = row
                added += 1
                log_entries.append({
                    'action': 'ADDED',
                    'gene_id': row.get('gene_id', ''),
                    'es_code': row.get('es_code', ''),
                    'source': filepath,
                    'note': '',
                })

    # Write merged output
    merged_rows = sorted(master_index.values(), key=lambda r: (r.get('es_code', ''), r.get('gene_id', '')))
    with open(args.output, 'w') as f:
        f.write('\t'.join(COLUMNS) + '\n')
        for row in merged_rows:
            f.write('\t'.join(str(row.get(c, '')) for c in COLUMNS) + '\n')

    # Write log
    log_cols = ['action', 'gene_id', 'es_code', 'source', 'note']
    with open(args.log, 'w') as f:
        f.write('\t'.join(log_cols) + '\n')
        for entry in log_entries:
            f.write('\t'.join(str(entry.get(c, '')) for c in log_cols) + '\n')

    print(f'Merge: {added} added, {skipped} skipped, {conflicts} role conflicts',
          file=sys.stderr)
    print(f'Total mappings: {len(merged_rows)}', file=sys.stderr)


if __name__ == '__main__':
    main()
