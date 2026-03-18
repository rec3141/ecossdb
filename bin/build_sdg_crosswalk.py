#!/usr/bin/env python3
"""Build the full CICES → MESH ES → SDG target crosswalk.

Chains the CICES-to-MESH-ES crosswalk with the MESH_SDG Lookup table
to produce a direct CICES code → SDG target mapping.

Usage:
    build_sdg_crosswalk.py \
        --cices-mesh db/ontology/sdg/cices_to_mesh_es.tsv \
        --mesh-sdg db/ontology/sdg/mesh_sdg_lookup.csv \
        --sdg-targets db/ontology/sdg/sdg_targets.tsv \
        --output db/ontology/sdg/cices_to_sdg.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict


def load_cices_mesh(filepath):
    """Load CICES → MESH ES mapping."""
    mapping = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            code = row['cices_code'].strip()
            mesh = row['mesh_es'].strip()
            if code not in mapping:
                mapping[code] = []
            mapping[code].append(mesh)
    return mapping


def load_mesh_sdg(filepath):
    """Load MESH ES → SDG Lookup (from MESH_SDG project).

    Returns dict: {mesh_es: [(goal, target, link_all, link_experts), ...]}
    """
    mapping = defaultdict(list)
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            es = row['ES'].strip()
            mapping[es].append({
                'sdg_goal': row['Goal_num'].strip(),
                'sdg_target': row['Target_num'].strip(),
                'link_all': row['LinkAll'].strip(),
                'link_experts': row['LinkExperts'].strip(),
            })
    return mapping


def load_sdg_targets(filepath):
    """Load SDG target descriptions."""
    targets = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            targets[row['sdg_target']] = {
                'goal_name': row['goal_name'],
                'target_description': row['target_description'],
            }
    return targets


def compute_link_strength(link_all, link_experts):
    """Convert MESH_SDG linkage codes to numeric strength.

    PO+PO = 1.0 (both agree positive)
    PO+UN or UN+PO = 0.5 (one source positive)
    UN+UN = 0.0 (no positive evidence, excluded)
    """
    if link_all == 'PO' and link_experts == 'PO':
        return 1.0
    elif link_all == 'PO' or link_experts == 'PO':
        return 0.5
    return 0.0


def main():
    parser = argparse.ArgumentParser(description='Build CICES → SDG crosswalk')
    parser.add_argument('--cices-mesh', required=True, help='CICES → MESH ES TSV')
    parser.add_argument('--mesh-sdg', required=True, help='MESH_SDG Lookup CSV')
    parser.add_argument('--sdg-targets', required=True, help='SDG targets TSV')
    parser.add_argument('--output', required=True, help='Output CICES → SDG TSV')
    parser.add_argument('--min-strength', type=float, default=0.5,
                        help='Minimum link strength to include (default: 0.5)')
    args = parser.parse_args()

    cices_mesh = load_cices_mesh(args.cices_mesh)
    mesh_sdg = load_mesh_sdg(args.mesh_sdg)
    sdg_targets = load_sdg_targets(args.sdg_targets)

    # Chain: CICES → MESH ES → SDG
    crosswalk = []
    for cices_code, mesh_labels in cices_mesh.items():
        for mesh_es in mesh_labels:
            if mesh_es not in mesh_sdg:
                continue
            for link in mesh_sdg[mesh_es]:
                strength = compute_link_strength(link['link_all'], link['link_experts'])
                if strength < args.min_strength:
                    continue

                target_info = sdg_targets.get(link['sdg_target'], {})
                crosswalk.append({
                    'cices_code': cices_code,
                    'mesh_es': mesh_es,
                    'sdg_goal': link['sdg_goal'],
                    'sdg_target': link['sdg_target'],
                    'goal_name': target_info.get('goal_name', ''),
                    'target_description': target_info.get('target_description', ''),
                    'link_strength': strength,
                    'link_all': link['link_all'],
                    'link_experts': link['link_experts'],
                })

    # Sort by CICES code then SDG target
    crosswalk.sort(key=lambda r: (r['cices_code'], r['sdg_target']))

    # Write output
    cols = ['cices_code', 'mesh_es', 'sdg_goal', 'sdg_target', 'goal_name',
            'target_description', 'link_strength', 'link_all', 'link_experts']
    with open(args.output, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in crosswalk:
            f.write('\t'.join(str(row[c]) for c in cols) + '\n')

    # Stats
    unique_cices = len(set(r['cices_code'] for r in crosswalk))
    unique_sdg = len(set(r['sdg_target'] for r in crosswalk))
    strong = sum(1 for r in crosswalk if r['link_strength'] == 1.0)
    print(f'Crosswalk: {len(crosswalk)} links ({strong} strong, '
          f'{len(crosswalk) - strong} moderate)', file=sys.stderr)
    print(f'Covering {unique_cices} CICES codes → {unique_sdg} SDG targets',
          file=sys.stderr)


if __name__ == '__main__':
    main()
