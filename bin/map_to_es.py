#!/usr/bin/env python3
"""Core mapping engine: gene annotations → ecosystem service assignments.

For each protein, checks annotation columns (KO, COG, EC, Pfam, CAZy)
against the gene-to-ES mapping table. One protein can map to multiple ES.

Usage:
    map_to_es.py \
        --annotations normalized_annotations.tsv \
        --mapping es_gene_mapping.tsv \
        --ontology es_ontology.tsv \
        --output es_gene_hits.tsv \
        --stats es_mapping_stats.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict


def load_mapping(mapping_path):
    """Load gene-to-ES mapping table.

    Returns dict: {(gene_id, gene_id_type): [mapping_row, ...]}
    """
    mapping = defaultdict(list)
    with open(mapping_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row['gene_id'], row['gene_id_type'])
            mapping[key].append(row)
    return mapping


def load_ontology(ontology_path):
    """Load ES ontology for name lookups."""
    ontology = {}
    with open(ontology_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ontology[row['es_code']] = row
    return ontology


# Map annotation columns to gene_id_type values
ANNOTATION_COLUMNS = [
    ('KO', 'KO'),
    ('COG', 'COG'),
    ('EC', 'EC'),
    ('Pfam', 'Pfam'),
    ('CAZy', 'CAZy'),
]


def map_protein(row, mapping):
    """Map a single protein to ES categories.

    Returns list of hit dicts.
    """
    hits = []
    protein_id = row['protein_id']
    contig_id = row.get('contig_id', '')

    for col, id_type in ANNOTATION_COLUMNS:
        value = row.get(col, '').strip()
        if not value:
            continue

        # Handle comma-separated multi-annotations
        gene_ids = [v.strip() for v in value.split(',') if v.strip()]

        for gene_id in gene_ids:
            key = (gene_id, id_type)
            if key in mapping:
                for m in mapping[key]:
                    hits.append({
                        'protein_id': protein_id,
                        'contig_id': contig_id,
                        'gene_id': gene_id,
                        'gene_id_type': id_type,
                        'es_code': m['es_code'],
                        'es_name': m['es_name'],
                        'confidence': m['confidence'],
                        'functional_role': m['functional_role'],
                        'evidence_source': m['evidence_source'],
                        'pathway_context': m.get('pathway_context', ''),
                        'detection_method': id_type.lower(),
                    })

    return hits


def main():
    parser = argparse.ArgumentParser(description='Map annotations to ES')
    parser.add_argument('--annotations', required=True, help='Normalized annotations TSV')
    parser.add_argument('--mapping', required=True, help='Gene-to-ES mapping table')
    parser.add_argument('--ontology', required=True, help='ES ontology TSV')
    parser.add_argument('--output', required=True, help='Output gene hits TSV')
    parser.add_argument('--stats', required=True, help='Output mapping stats TSV')
    args = parser.parse_args()

    mapping = load_mapping(args.mapping)
    ontology = load_ontology(args.ontology)

    # Process annotations
    total_proteins = 0
    mapped_proteins = 0
    all_hits = []

    with open(args.annotations) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            total_proteins += 1
            hits = map_protein(row, mapping)
            if hits:
                mapped_proteins += 1
                all_hits.extend(hits)

    # Write hits
    hit_cols = ['protein_id', 'contig_id', 'gene_id', 'gene_id_type',
                'es_code', 'es_name', 'confidence', 'functional_role',
                'evidence_source', 'pathway_context', 'detection_method']

    with open(args.output, 'w') as f:
        f.write('\t'.join(hit_cols) + '\n')
        for hit in all_hits:
            f.write('\t'.join(str(hit.get(c, '')) for c in hit_cols) + '\n')

    # Write stats
    es_counts = defaultdict(int)
    for hit in all_hits:
        es_counts[hit['es_code']] += 1

    with open(args.stats, 'w') as f:
        f.write('\t'.join(['metric', 'value']) + '\n')
        f.write(f'total_proteins\t{total_proteins}\n')
        f.write(f'mapped_proteins\t{mapped_proteins}\n')
        f.write(f'total_hits\t{len(all_hits)}\n')
        f.write(f'unique_es_codes\t{len(es_counts)}\n')
        pct = (mapped_proteins / total_proteins * 100) if total_proteins else 0
        f.write(f'mapping_rate_pct\t{pct:.1f}\n')
        f.write('\n')
        f.write('\t'.join(['es_code', 'es_name', 'hit_count']) + '\n')
        for code in sorted(es_counts, key=lambda c: -es_counts[c]):
            name = ontology.get(code, {}).get('name', '')
            f.write(f'{code}\t{name}\t{es_counts[code]}\n')

    print(f'Mapped {mapped_proteins}/{total_proteins} proteins '
          f'({pct:.1f}%) → {len(all_hits)} ES hits across {len(es_counts)} categories',
          file=sys.stderr)


if __name__ == '__main__':
    main()
