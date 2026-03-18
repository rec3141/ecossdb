#!/usr/bin/env python3
"""Parse hmmsearch domtblout output and map hits to ES categories.

Applies per-profile thresholds from metadata TSV and maps hits
to ecosystem services using HMM metadata.

Usage:
    hmm_to_es.py \
        --domtblout hmmsearch_results.domtblout \
        --metadata es_hmm_metadata.tsv \
        --output hmm_es_hits.tsv
"""

import argparse
import csv
import sys


def load_metadata(filepath):
    """Load HMM metadata: profile_name → {es_code, es_name, confidence, ...}"""
    meta = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            meta[row['profile_name']] = {
                'es_code': row['es_code'],
                'es_name': row['es_name'],
                'confidence': float(row.get('confidence', 0.7)),
                'functional_role': row.get('functional_role', 'producer'),
                'score_threshold': float(row.get('score_threshold', 0)),
                'evalue_threshold': float(row.get('evalue_threshold', 1e-5)),
                'description': row.get('description', ''),
            }
    return meta


def parse_domtblout(filepath, metadata):
    """Parse hmmsearch domtblout and filter by per-profile thresholds."""
    hits = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue

            protein_id = parts[0]
            profile_name = parts[3]
            evalue = float(parts[6])
            score = float(parts[7])

            if profile_name not in metadata:
                continue

            meta = metadata[profile_name]
            if score < meta['score_threshold']:
                continue
            if evalue > meta['evalue_threshold']:
                continue

            contig_id = '_'.join(protein_id.rsplit('_', 1)[:-1]) if '_' in protein_id else protein_id

            hits.append({
                'protein_id': protein_id,
                'contig_id': contig_id,
                'gene_id': profile_name,
                'gene_id_type': 'custom_hmm',
                'es_code': meta['es_code'],
                'es_name': meta['es_name'],
                'confidence': meta['confidence'],
                'functional_role': meta['functional_role'],
                'evidence_source': 'custom_hmm',
                'pathway_context': meta['description'],
                'detection_method': 'hmmsearch',
                'hmm_score': score,
                'hmm_evalue': evalue,
            })

    return hits


def main():
    parser = argparse.ArgumentParser(description='Map HMM hits to ES')
    parser.add_argument('--domtblout', required=True, help='hmmsearch domtblout file')
    parser.add_argument('--metadata', required=True, help='HMM metadata TSV')
    parser.add_argument('--output', required=True, help='Output ES hits TSV')
    args = parser.parse_args()

    metadata = load_metadata(args.metadata)
    hits = parse_domtblout(args.domtblout, metadata)

    cols = ['protein_id', 'contig_id', 'gene_id', 'gene_id_type',
            'es_code', 'es_name', 'confidence', 'functional_role',
            'evidence_source', 'pathway_context', 'detection_method']
    with open(args.output, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for hit in hits:
            f.write('\t'.join(str(hit.get(c, '')) for c in cols) + '\n')

    print(f'Mapped {len(hits)} HMM hits from {len(metadata)} profiles', file=sys.stderr)


if __name__ == '__main__':
    main()
