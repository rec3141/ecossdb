#!/usr/bin/env python3
"""Aggregate ES gene hits to contig level and produce gene catalog.

Merges annotation-based and optional HMM-based hits, deduplicates
(higher confidence wins for same protein+ES), and rolls up to contigs.

Usage:
    aggregate_contigs.py \
        --gene-hits es_gene_hits.tsv \
        [--hmm-hits hmm_es_hits.tsv] \
        --output-catalog es_gene_catalog.tsv \
        --output-contigs es_per_contig.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict


def load_hits(filepath):
    """Load gene-ES hits from TSV."""
    hits = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            row['confidence'] = float(row.get('confidence', 0))
            hits.append(row)
    return hits


def merge_and_deduplicate(annotation_hits, hmm_hits):
    """Merge annotation and HMM hits, keeping highest confidence per protein+ES."""
    # Key: (protein_id, es_code)
    best = {}
    for hit in annotation_hits + hmm_hits:
        key = (hit['protein_id'], hit['es_code'])
        if key not in best or hit['confidence'] > best[key]['confidence']:
            best[key] = hit
    return sorted(best.values(), key=lambda h: (h.get('contig_id', ''), h['protein_id']))


def aggregate_to_contigs(catalog):
    """Roll up gene catalog to per-contig ES profiles."""
    contigs = defaultdict(lambda: defaultdict(lambda: {
        'gene_count': 0,
        'max_confidence': 0,
        'sum_confidence': 0,
        'roles': defaultdict(int),
        'genes': [],
    }))

    for hit in catalog:
        contig_id = hit.get('contig_id', '')
        if not contig_id:
            continue
        es_code = hit['es_code']
        entry = contigs[contig_id][es_code]
        entry['gene_count'] += 1
        conf = hit['confidence']
        entry['max_confidence'] = max(entry['max_confidence'], conf)
        entry['sum_confidence'] += conf
        entry['roles'][hit.get('functional_role', 'unknown')] += 1
        entry['genes'].append(f"{hit.get('gene_id', '')}({conf:.2f})")
        entry['es_name'] = hit.get('es_name', '')

    return contigs


def main():
    parser = argparse.ArgumentParser(description='Aggregate ES hits to contigs')
    parser.add_argument('--gene-hits', required=True, help='Annotation-based ES hits')
    parser.add_argument('--hmm-hits', help='HMM-based ES hits (optional)')
    parser.add_argument('--output-catalog', required=True, help='Output gene catalog TSV')
    parser.add_argument('--output-contigs', required=True, help='Output per-contig TSV')
    args = parser.parse_args()

    annotation_hits = load_hits(args.gene_hits)
    hmm_hits = load_hits(args.hmm_hits) if args.hmm_hits else []

    catalog = merge_and_deduplicate(annotation_hits, hmm_hits)

    # Write gene catalog
    cat_cols = ['protein_id', 'contig_id', 'gene_id', 'gene_id_type',
                'es_code', 'es_name', 'confidence', 'functional_role',
                'detection_method']
    with open(args.output_catalog, 'w') as f:
        f.write('\t'.join(cat_cols) + '\n')
        for hit in catalog:
            f.write('\t'.join(str(hit.get(c, '')) for c in cat_cols) + '\n')

    # Aggregate and write per-contig
    contigs = aggregate_to_contigs(catalog)
    contig_cols = ['contig_id', 'es_code', 'es_name', 'gene_count',
                   'mean_confidence', 'max_confidence', 'functional_roles', 'top_genes']

    with open(args.output_contigs, 'w') as f:
        f.write('\t'.join(contig_cols) + '\n')
        for contig_id in sorted(contigs):
            for es_code in sorted(contigs[contig_id]):
                entry = contigs[contig_id][es_code]
                mean_conf = entry['sum_confidence'] / entry['gene_count']
                roles_str = ','.join(f'{r}:{c}' for r, c in sorted(entry['roles'].items()))
                genes_str = ','.join(entry['genes'][:10])  # top 10
                f.write('\t'.join([
                    contig_id, es_code, entry['es_name'],
                    str(entry['gene_count']),
                    f'{mean_conf:.3f}', f'{entry["max_confidence"]:.3f}',
                    roles_str, genes_str,
                ]) + '\n')

    print(f'Gene catalog: {len(catalog)} entries across {len(contigs)} contigs',
          file=sys.stderr)


if __name__ == '__main__':
    main()
