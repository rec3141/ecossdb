#!/usr/bin/env python3
"""Generate visualization JSON for the danaseq Ecosystem Services tab.

Produces ecosystem_services.json (+ .gz) containing all data needed
for treemap, sankey, heatmap, and detail table panels.

Usage:
    es_to_json.py \
        --scores es_scores.tsv \
        --catalog es_gene_catalog.tsv \
        --mag-profiles es_per_mag.tsv \
        --hierarchy es_hierarchy.json \
        --output ecosystem_services.json
"""

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict


def load_tsv(filepath):
    """Load TSV into list of dicts."""
    with open(filepath) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def build_treemap_data(scores, hierarchy):
    """Build treemap data: CICES hierarchy sized by weighted score."""
    # Aggregate scores across all entities
    es_totals = defaultdict(float)
    es_confidence = defaultdict(list)
    for row in scores:
        code = row['es_code']
        es_totals[code] += float(row.get('weighted_score', 0))
        es_confidence[code].append(float(row.get('completeness', 0)))

    def annotate_tree(node):
        code = node.get('code', '')
        if code in es_totals:
            node['value'] = round(es_totals[code], 4)
            confs = es_confidence[code]
            node['confidence'] = round(sum(confs) / len(confs), 3) if confs else 0
        if 'children' in node:
            for child in node['children']:
                annotate_tree(child)
            # Propagate values up
            child_sum = sum(c.get('value', 0) for c in node['children'])
            if child_sum > 0 and 'value' not in node:
                node['value'] = round(child_sum, 4)

    annotate_tree(hierarchy)
    return hierarchy


def build_sankey_data(catalog):
    """Build sankey data: functional categories → ES categories."""
    links = defaultdict(int)
    for row in catalog:
        pathway = row.get('pathway_context', row.get('detection_method', 'unknown'))
        # Simplify pathway to category
        category = pathway.split('/')[0].split('(')[0].strip()
        if not category:
            category = 'Other'
        es_name = row.get('es_name', row.get('es_code', ''))
        role = row.get('functional_role', 'unknown')
        links[(category, es_name, role)] += 1

    nodes = set()
    sankey_links = []
    for (source, target, role), count in links.items():
        nodes.add(source)
        nodes.add(target)
        sankey_links.append({
            'source': source,
            'target': target,
            'value': count,
            'role': role,
        })

    return {
        'nodes': sorted(nodes),
        'links': sankey_links,
    }


def build_heatmap_data(mag_profiles):
    """Build heatmap data: MAGs x ES codes."""
    mags = set()
    es_codes = set()
    values = {}

    for row in mag_profiles:
        mag = row['entity_id']
        es = row['es_code']
        mags.add(mag)
        es_codes.add(es)
        values[(mag, es)] = float(row.get('weighted_score', 0))

    mag_list = sorted(mags)
    es_list = sorted(es_codes)

    matrix = []
    for mag in mag_list:
        row_vals = [round(values.get((mag, es), 0), 4) for es in es_list]
        matrix.append(row_vals)

    return {
        'mags': mag_list,
        'es_codes': es_list,
        'es_names': {row['es_code']: row.get('es_name', '') for row in mag_profiles},
        'matrix': matrix,
    }


def main():
    parser = argparse.ArgumentParser(description='Generate viz JSON')
    parser.add_argument('--scores', required=True, help='ES scores TSV')
    parser.add_argument('--catalog', required=True, help='Gene catalog TSV')
    parser.add_argument('--mag-profiles', required=True, help='Per-MAG profiles TSV')
    parser.add_argument('--hierarchy', required=True, help='ES hierarchy JSON')
    parser.add_argument('--output', required=True, help='Output JSON path')
    args = parser.parse_args()

    scores = load_tsv(args.scores)
    catalog = load_tsv(args.catalog)
    mag_profiles = load_tsv(args.mag_profiles)

    with open(args.hierarchy) as f:
        hierarchy = json.load(f)

    output = {
        'treemap': build_treemap_data(scores, hierarchy),
        'sankey': build_sankey_data(catalog),
        'heatmap': build_heatmap_data(mag_profiles),
        'catalog': catalog[:1000],  # cap detail table at 1000 rows
    }

    with open(args.output, 'w') as f:
        json.dump(output, f)

    # Also write gzipped version
    gz_path = args.output + '.gz'
    with gzip.open(gz_path, 'wt') as f:
        json.dump(output, f)

    print(f'Wrote viz JSON: {len(scores)} scores, {len(catalog)} catalog entries, '
          f'{len(mag_profiles)} MAG profiles', file=sys.stderr)


if __name__ == '__main__':
    main()
