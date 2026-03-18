#!/usr/bin/env python3
"""Aggregate MAG-level ES profiles to cross-sample ES matrix.

Produces a samples x ES_codes matrix with weighted scores.

Usage:
    aggregate_samples.py \
        --mag-profiles es_per_mag_*.tsv \
        --output-matrix es_matrix.tsv \
        --output-summary es_summary.tsv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='Cross-sample ES aggregation')
    parser.add_argument('--mag-profiles', required=True, nargs='+',
                        help='Per-MAG ES profile TSVs (one per sample)')
    parser.add_argument('--sample-ids', nargs='+',
                        help='Sample IDs (default: derived from filenames)')
    parser.add_argument('--output-matrix', required=True, help='Output ES matrix')
    parser.add_argument('--output-summary', required=True, help='Output ES summary')
    args = parser.parse_args()

    # Determine sample IDs
    if args.sample_ids and len(args.sample_ids) == len(args.mag_profiles):
        sample_ids = args.sample_ids
    else:
        sample_ids = [os.path.basename(f).replace('es_per_mag_', '').replace('.tsv', '')
                      for f in args.mag_profiles]

    # Aggregate: sample → es_code → score
    matrix = defaultdict(lambda: defaultdict(float))
    all_es_codes = set()
    es_names = {}

    for sample_id, filepath in zip(sample_ids, args.mag_profiles):
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                es_code = row['es_code']
                score = float(row.get('weighted_score', 0))
                matrix[sample_id][es_code] += score
                all_es_codes.add(es_code)
                es_names[es_code] = row.get('es_name', '')

    es_codes_sorted = sorted(all_es_codes)

    # Write matrix
    with open(args.output_matrix, 'w') as f:
        f.write('\t'.join(['sample_id'] + es_codes_sorted) + '\n')
        for sample_id in sample_ids:
            values = [f'{matrix[sample_id].get(ec, 0):.4f}' for ec in es_codes_sorted]
            f.write('\t'.join([sample_id] + values) + '\n')

    # Write summary
    with open(args.output_summary, 'w') as f:
        f.write('\t'.join(['es_code', 'es_name', 'n_samples', 'mean_score',
                           'max_score', 'prevalence']) + '\n')
        for ec in es_codes_sorted:
            scores = [matrix[s][ec] for s in sample_ids if ec in matrix[s]]
            n = len(scores)
            mean_s = sum(scores) / n if n else 0
            max_s = max(scores) if scores else 0
            prev = n / len(sample_ids) if sample_ids else 0
            f.write('\t'.join([
                ec, es_names.get(ec, ''), str(n),
                f'{mean_s:.4f}', f'{max_s:.4f}', f'{prev:.3f}',
            ]) + '\n')

    print(f'Matrix: {len(sample_ids)} samples x {len(es_codes_sorted)} ES codes',
          file=sys.stderr)


if __name__ == '__main__':
    main()
