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

DEFAULT_ROLE_WEIGHTS = {
    'producer': 1.0,
    'transformer': 0.7,
    'consumer': -0.3,
    'inhibitor': -0.5,
}


def load_contig2bin(filepath):
    """Load contig-to-bin mapping. Format: contig_id<TAB>bin_id"""
    mapping = {}
    with open(filepath) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 2 and not row[0].startswith('#'):
                mapping[row[0].strip()] = row[1].strip()
    return mapping


def parse_role_weights(weights_str):
    """Parse role:weight pairs from comma-separated string."""
    weights = {}
    for pair in weights_str.split(','):
        role, weight = pair.strip().split(':')
        weights[role.strip()] = float(weight.strip())
    return weights


def main():
    parser = argparse.ArgumentParser(description='Aggregate ES profiles per MAG')
    parser.add_argument('--catalog', required=True, help='Gene catalog TSV')
    parser.add_argument('--contig2bin', required=True, help='Contig-to-bin mapping')
    parser.add_argument('--output', required=True, help='Output per-MAG ES profile')
    parser.add_argument('--role-weights',
                        default='producer:1.0,transformer:0.7,consumer:-0.3,inhibitor:-0.5',
                        help='Role weights for scoring')
    args = parser.parse_args()

    contig2bin = load_contig2bin(args.contig2bin)
    role_weights = parse_role_weights(args.role_weights)

    # Aggregate catalog entries per bin
    bins = defaultdict(lambda: defaultdict(lambda: {
        'gene_count': 0,
        'sum_weighted': 0,
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
            role = row.get('functional_role', 'unknown')
            weight = role_weights.get(role, 0)
            entry['sum_weighted'] += conf * weight
            entry['sum_confidence'] += conf
            entry['max_confidence'] = max(entry['max_confidence'], conf)
            entry['roles'][role] += 1
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
                # Role-weighted mean: mean(confidence × role_weight)
                mean_weighted = entry['sum_weighted'] / entry['gene_count']
                roles_str = ','.join(f'{r}:{c}' for r, c in sorted(entry['roles'].items()))
                # Top 5 genes by confidence
                top = sorted(entry['genes'], key=lambda g: -float(g.split('(')[1].rstrip(')')))[:5]
                f.write('\t'.join([
                    bin_id, 'MAG', es_code, entry['es_name'],
                    str(entry['gene_count']),
                    f'{mean_weighted:.4f}',
                    '-',  # completeness filled by score_es.py
                    roles_str,
                    ','.join(top),
                ]) + '\n')

    binned = total - unbinned
    print(f'Aggregated {binned}/{total} genes into {len(bins)} MAGs '
          f'({unbinned} unbinned)', file=sys.stderr)


if __name__ == '__main__':
    main()
