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
import math
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


def cluster_order(matrix):
    """Hierarchical clustering using Bray-Curtis distance.

    Uses scipy for fast computation (handles 8K+ rows).
    Falls back to greedy nearest-neighbour in pure Python if scipy unavailable.
    Returns a permutation of row indices.
    """
    n = len(matrix)
    if n <= 2:
        return list(range(n))

    try:
        import numpy as np
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, leaves_list

        arr = np.array(matrix, dtype=np.float64)
        # Remove all-zero rows for distance computation (BC undefined)
        row_sums = arr.sum(axis=1)
        nonzero = row_sums > 0
        if nonzero.sum() < 2:
            return list(range(n))

        # Compute on nonzero rows only
        nz_indices = np.where(nonzero)[0]
        nz_arr = arr[nz_indices]
        dists = pdist(nz_arr, metric='braycurtis')
        # Replace NaN with 1.0 (max dissimilarity)
        dists = np.nan_to_num(dists, nan=1.0)
        Z = linkage(dists, method='average')
        nz_order = leaves_list(Z)
        # Map back to original indices, append zero rows at the end
        zero_indices = np.where(~nonzero)[0]
        order = [int(nz_indices[i]) for i in nz_order] + [int(i) for i in zero_indices]
        return order

    except ImportError:
        # Pure Python fallback: greedy nearest-neighbour
        def bray_curtis(a, b):
            num = sum(abs(x - y) for x, y in zip(a, b))
            den = sum(x + y for x, y in zip(a, b))
            return num / den if den > 0 else 0.0

        dist = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d = bray_curtis(matrix[i], matrix[j])
                dist[i][j] = d
                dist[j][i] = d

        sums = [sum(row) for row in matrix]
        start = max(range(n), key=lambda i: sums[i])
        visited = [False] * n
        order = [start]
        visited[start] = True
        for _ in range(n - 1):
            curr = order[-1]
            best_j, best_d = -1, float('inf')
            for j in range(n):
                if not visited[j] and dist[curr][j] < best_d:
                    best_d = dist[curr][j]
                    best_j = j
            order.append(best_j)
            visited[best_j] = True
        return order


def build_pathway_heatmap(catalog, mapping_lookup, contig2bin):
    """Build expanded heatmap: bins x (ES code + pathway_context).

    contig2bin may map contigs to multiple bins (from different binners).
    One contig can appear in bins from each binner.
    """
    # Count hits per bin per (es_code, pathway)
    bin_pathway = defaultdict(lambda: defaultdict(int))
    for row in catalog:
        contig = row.get('contig_id', '')
        bins = contig2bin.get(contig, [])
        if not bins:
            continue
        gene_id = row['gene_id']
        es_code = row['es_code']
        pw = mapping_lookup.get((gene_id, es_code), {}).get('pathway_context', '')
        if not pw:
            pw = 'Other'
        key = (es_code, pw)
        for mag in bins:
            bin_pathway[mag][key] += 1

    # Build sorted column list
    all_cols = set()
    for mag_data in bin_pathway.values():
        all_cols.update(mag_data.keys())
    col_list = sorted(all_cols)

    mag_list = sorted(bin_pathway.keys())
    matrix = []
    for mag in mag_list:
        row_vals = [bin_pathway[mag].get(col, 0) for col in col_list]
        matrix.append(row_vals)

    # Bray-Curtis column clustering
    col_order = list(range(len(col_list)))
    if matrix and matrix[0]:
        transposed = [[matrix[r][c] for r in range(len(matrix))] for c in range(len(col_list))]
        print(f'  Computing Bray-Curtis column clustering ({len(col_list)} cols)...', file=sys.stderr)
        col_order = cluster_order(transposed)

    # Full Bray-Curtis row clustering (scipy handles 8K+ rows efficiently)
    print(f'  Computing Bray-Curtis row clustering ({len(mag_list)} bins)...', file=sys.stderr)
    row_order = cluster_order(matrix)

    # Extract binner name per bin
    def get_binner(name):
        if name.startswith('dastool-'):
            return 'dastool'
        return name.rsplit('_', 1)[0] if '_' in name else name

    return {
        'mags': mag_list,
        'columns': [{'es_code': c[0], 'pathway': c[1]} for c in col_list],
        'matrix': matrix,
        'row_order': row_order,
        'col_order': col_order,
        'binners': [get_binner(m) for m in mag_list],
    }


def main():
    parser = argparse.ArgumentParser(description='Generate viz JSON')
    parser.add_argument('--scores', required=True, help='ES scores TSV')
    parser.add_argument('--catalog', required=True, help='Gene catalog TSV')
    parser.add_argument('--mag-profiles', required=True, help='Per-MAG profiles TSV')
    parser.add_argument('--contig2bin', default=None, help='Contig-to-bin mapping TSV')
    parser.add_argument('--hierarchy', required=True, help='ES hierarchy JSON')
    parser.add_argument('--mapping', default=None, help='ES gene mapping TSV (for pathway_context)')
    parser.add_argument('--output', required=True, help='Output JSON path')
    args = parser.parse_args()

    scores = load_tsv(args.scores)
    catalog = load_tsv(args.catalog)
    mag_profiles = load_tsv(args.mag_profiles)

    # Load mapping table for pathway_context enrichment
    mapping_lookup = {}
    if args.mapping:
        for row in load_tsv(args.mapping):
            mapping_lookup[(row['gene_id'], row['es_code'])] = {
                'pathway_context': row.get('pathway_context', ''),
                'notes': row.get('notes', ''),
            }

    with open(args.hierarchy) as f:
        hierarchy = json.load(f)

    # Aggregate catalog to unique gene × ES entries with hit counts
    gene_es = {}
    for row in catalog:
        key = (row['gene_id'], row['es_code'])
        if key not in gene_es:
            gene_es[key] = {
                'gene_id': row['gene_id'],
                'gene_id_type': row.get('gene_id_type', ''),
                'es_code': row['es_code'],
                'es_name': row.get('es_name', ''),
                'confidence': row.get('confidence', ''),
                'functional_role': row.get('functional_role', ''),
                'detection_method': row.get('detection_method', ''),
                'pathway_context': mapping_lookup.get(key, {}).get('pathway_context', ''),
                'notes': mapping_lookup.get(key, {}).get('notes', ''),
                'n_hits': 0,
            }
        gene_es[key]['n_hits'] += 1
    # Round confidence to avoid floating point noise
    for entry in gene_es.values():
        try:
            entry['confidence'] = round(float(entry['confidence']), 2)
        except (ValueError, TypeError):
            pass
    gene_catalog = sorted(gene_es.values(), key=lambda r: (-r['n_hits'], r['es_code'], r['gene_id']))

    # Build pathway-level heatmap if contig2bin provided
    # contig2bin maps contig → list of bins (one contig can be in multiple binner bins)
    contig2bin = defaultdict(list)
    if args.contig2bin:
        for c2b_path in args.contig2bin.split(','):
            c2b_path = c2b_path.strip()
            if not c2b_path:
                continue
            with open(c2b_path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2 and not parts[0].startswith('#'):
                        contig2bin[parts[0]].append(parts[1])

    pathway_heatmap = None
    if contig2bin and mapping_lookup:
        pathway_heatmap = build_pathway_heatmap(catalog, mapping_lookup, contig2bin)

    output = {
        'treemap': build_treemap_data(scores, hierarchy),
        'sankey': build_sankey_data(catalog),
        'heatmap': build_heatmap_data(mag_profiles),
        'catalog': gene_catalog,
        'total_protein_hits': len(catalog),
    }
    if pathway_heatmap:
        output['pathway_heatmap'] = pathway_heatmap

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
