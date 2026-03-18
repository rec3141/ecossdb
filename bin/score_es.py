#!/usr/bin/env python3
"""Confidence-weighted ES scoring with pathway completeness.

Composite scoring: gene_presence x confidence x role_weight x pathway_completeness

Role weights (default):
  producer:    +1.0
  transformer: +0.7
  consumer:    -0.3
  inhibitor:   -0.5

Usage:
    score_es.py \
        --catalog es_gene_catalog.tsv \
        --mapping es_gene_mapping.tsv \
        --role-weights producer:1.0,transformer:0.7,consumer:-0.3,inhibitor:-0.5 \
        --output es_scores.tsv \
        --confidence-out es_confidence.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict


def parse_role_weights(weights_str):
    """Parse role:weight pairs from comma-separated string."""
    weights = {}
    for pair in weights_str.split(','):
        role, weight = pair.strip().split(':')
        weights[role.strip()] = float(weight.strip())
    return weights


def load_pathway_steps(mapping_path):
    """Build pathway step sets from mapping table for completeness scoring.

    Groups genes by (es_code, pathway_context) to define "steps" analogous
    to KEGG module steps. Each unique pathway_context within an ES code
    is one step; genes within that step are alternatives.
    """
    pathways = defaultdict(lambda: defaultdict(set))
    with open(mapping_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            es_code = row['es_code']
            pathway = row.get('pathway_context', '').strip()
            if pathway:
                pathways[es_code][pathway].add(row['gene_id'])
    return pathways


def evaluate_completeness(gene_set, pathway_steps):
    """Return fraction of pathway steps satisfied (0.0 to 1.0)."""
    if not pathway_steps:
        return 1.0  # no pathway info = assume complete
    present = sum(1 for step_genes in pathway_steps.values() if gene_set & step_genes)
    return present / len(pathway_steps)


def main():
    parser = argparse.ArgumentParser(description='Score ecosystem services')
    parser.add_argument('--catalog', required=True, help='Gene catalog TSV')
    parser.add_argument('--mapping', required=True, help='Gene-to-ES mapping table')
    parser.add_argument('--role-weights', default='producer:1.0,transformer:0.7,consumer:-0.3,inhibitor:-0.5')
    parser.add_argument('--output', required=True, help='Output ES scores TSV')
    parser.add_argument('--confidence-out', required=True, help='Output confidence detail TSV')
    parser.add_argument('--level', default='bin', choices=['contig', 'bin', 'sample'],
                        help='Aggregation level (default: bin)')
    args = parser.parse_args()

    role_weights = parse_role_weights(args.role_weights)
    pathway_steps = load_pathway_steps(args.mapping)

    # Determine entity column
    entity_col = {
        'contig': 'contig_id',
        'bin': 'bin_id',
        'sample': 'sample_id',
    }.get(args.level, 'contig_id')

    # Aggregate genes per entity per ES code
    entities = defaultdict(lambda: defaultdict(lambda: {
        'genes': set(),
        'gene_details': [],
        'roles': defaultdict(int),
        'sum_weighted': 0,
        'es_name': '',
    }))

    with open(args.catalog) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            entity_id = row.get(entity_col, row.get('contig_id', ''))
            if not entity_id:
                continue
            es_code = row['es_code']
            entry = entities[entity_id][es_code]
            gene_id = row.get('gene_id', '')
            conf = float(row.get('confidence', 0))
            role = row.get('functional_role', 'unknown')
            weight = role_weights.get(role, 0)

            entry['genes'].add(gene_id)
            entry['gene_details'].append((gene_id, conf, role))
            entry['roles'][role] += 1
            entry['sum_weighted'] += conf * weight
            entry['es_name'] = row.get('es_name', '')

    # Score each entity x ES
    score_rows = []
    confidence_rows = []

    for entity_id in sorted(entities):
        for es_code in sorted(entities[entity_id]):
            entry = entities[entity_id][es_code]
            gene_set = entry['genes']
            steps = pathway_steps.get(es_code, {})
            completeness = evaluate_completeness(gene_set, steps)

            n_genes = len(gene_set)
            raw_score = entry['sum_weighted'] / n_genes if n_genes else 0
            final_score = raw_score * completeness

            roles_str = ','.join(f'{r}:{c}' for r, c in sorted(entry['roles'].items()))

            score_rows.append({
                'entity_id': entity_id,
                'es_code': es_code,
                'es_name': entry['es_name'],
                'gene_count': n_genes,
                'weighted_score': f'{final_score:.4f}',
                'raw_score': f'{raw_score:.4f}',
                'completeness': f'{completeness:.3f}',
                'functional_roles': roles_str,
            })

            for gene_id, conf, role in entry['gene_details']:
                confidence_rows.append({
                    'entity_id': entity_id,
                    'es_code': es_code,
                    'gene_id': gene_id,
                    'confidence': f'{conf:.3f}',
                    'functional_role': role,
                    'role_weight': f'{role_weights.get(role, 0):.2f}',
                })

    # Write scores
    score_cols = ['entity_id', 'es_code', 'es_name', 'gene_count',
                  'weighted_score', 'raw_score', 'completeness', 'functional_roles']
    with open(args.output, 'w') as f:
        f.write('\t'.join(score_cols) + '\n')
        for row in score_rows:
            f.write('\t'.join(str(row[c]) for c in score_cols) + '\n')

    # Write confidence detail
    conf_cols = ['entity_id', 'es_code', 'gene_id', 'confidence',
                 'functional_role', 'role_weight']
    with open(args.confidence_out, 'w') as f:
        f.write('\t'.join(conf_cols) + '\n')
        for row in confidence_rows:
            f.write('\t'.join(str(row[c]) for c in conf_cols) + '\n')

    print(f'Scored {len(score_rows)} entity-ES pairs across {len(entities)} entities',
          file=sys.stderr)


if __name__ == '__main__':
    main()
