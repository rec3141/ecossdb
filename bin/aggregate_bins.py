#!/usr/bin/env python3
"""Aggregate ES profiles from contig level to MAG/bin level.

Uses contig-to-bin mapping (same format as DAS Tool output) to roll up
per-contig ES hits into per-MAG profiles.

Usage:
    aggregate_bins.py \
        --catalog es_gene_catalog.tsv \
        --contig2bin contig2bin.tsv \
        --output es_per_mag.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict


def load_contig2bin(filepath):
    """Load contig-to-bin mapping. Format: contig_id<TAB>bin_id"""
    mapping = {}
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 2 and not row[0].startswith('#'):
                mapping[row[0].strip()] = row[1].strip()
    return mapping


def main():
    parser = argparse.ArgumentParser(description='Aggregate ES profiles per MAG')
    parser.add_argument('--catalog', required=True, help='Gene catalog TSV')
    parser.add_argument('--contig2bin', required=True, help='Contig-to-bin mapping')
    parser.add_argument('--output', required=True, help='Output per-MAG ES profile')
    args = parser.parse_args()

    contig2bin = load_contig2bin(args.contig2bin)

    # Aggregate catalog entries per bin
    bins = defaultdict(lambda: defaultdict(lambda: {
        'gene_count': 0,
        'sum_confidence': 0,
        'max_confidence': 0,
        'roles': defaultdict(int),
        'genes': [],
        'es_name': '',
    }))

    unbinned = 0
    total = 0

    with open(args.catalog) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total += 1
            contig_id = row.get('contig_id', '')
            bin_id = contig2bin.get(contig_id)
            if not bin_id:
                unbinned += 1
                continue

            es_code = row['es_code']
            entry = bins[bin_id][es_code]
            entry['gene_count'] += 1
            conf = float(row.get('confidence', 0))
            entry['sum_confidence'] += conf
            entry['max_confidence'] = max(entry['max_confidence'], conf)
            entry['roles'][row.get('functional_role', 'unknown')] += 1
            entry['genes'].append(f"{row.get('gene_id', '')}({conf:.2f})")
            entry['es_name'] = row.get('es_name', '')

    # Write per-MAG profiles
    out_cols = ['entity_id', 'entity_type', 'es_code', 'es_name', 'gene_count',
                'weighted_score', 'completeness', 'functional_roles', 'top_genes']

    with open(args.output, 'w') as f:
        f.write('\t'.join(out_cols) + '\n')
        for bin_id in sorted(bins):
            for es_code in sorted(bins[bin_id]):
                entry = bins[bin_id][es_code]
                mean_conf = entry['sum_confidence'] / entry['gene_count']
                roles_str = ','.join(f'{r}:{c}' for r, c in sorted(entry['roles'].items()))
                # Top 5 genes by confidence
                top = sorted(entry['genes'], key=lambda g: -float(g.split('(')[1].rstrip(')')))[:5]
                f.write('\t'.join([
                    bin_id, 'MAG', es_code, entry['es_name'],
                    str(entry['gene_count']),
                    f'{mean_conf:.3f}',
                    '-',  # completeness filled by score_es.py
                    roles_str,
                    ','.join(top),
                ]) + '\n')

    binned = total - unbinned
    print(f'Aggregated {binned}/{total} genes into {len(bins)} MAGs '
          f'({unbinned} unbinned)', file=sys.stderr)


if __name__ == '__main__':
    main()
