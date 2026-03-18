#!/usr/bin/env python3
"""Map ES scores to SDG target contributions.

Takes ES scores (from score_es.py) and the CICES→SDG crosswalk to produce
SDG-level aggregations showing which Sustainable Development Goals a
microbial community supports.

Usage:
    map_es_to_sdg.py \
        --scores es_scores.tsv \
        --crosswalk db/ontology/sdg/cices_to_sdg.tsv \
        --sdg-targets db/ontology/sdg/sdg_targets.tsv \
        --output-targets es_sdg_targets.tsv \
        --output-goals es_sdg_goals.tsv \
        --output-json es_sdg.json
"""

import argparse
import csv
import json
import sys
from collections import defaultdict


def load_crosswalk(filepath):
    """Load CICES → SDG crosswalk."""
    crosswalk = defaultdict(list)  # cices_code → [sdg links]
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            crosswalk[row['cices_code']].append({
                'sdg_goal': row['sdg_goal'],
                'sdg_target': row['sdg_target'],
                'goal_name': row['goal_name'],
                'target_description': row.get('target_description', ''),
                'link_strength': float(row['link_strength']),
                'mesh_es': row['mesh_es'],
            })
    return crosswalk


def load_sdg_targets(filepath):
    """Load SDG target metadata."""
    targets = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            targets[row['sdg_target']] = row
    return targets


def main():
    parser = argparse.ArgumentParser(description='Map ES scores to SDG targets')
    parser.add_argument('--scores', required=True, help='ES scores TSV')
    parser.add_argument('--crosswalk', required=True, help='CICES → SDG crosswalk TSV')
    parser.add_argument('--sdg-targets', required=True, help='SDG targets TSV')
    parser.add_argument('--output-targets', required=True, help='Output per-target SDG scores')
    parser.add_argument('--output-goals', required=True, help='Output per-goal SDG scores')
    parser.add_argument('--output-json', required=True, help='Output JSON for viz')
    args = parser.parse_args()

    crosswalk = load_crosswalk(args.crosswalk)
    sdg_meta = load_sdg_targets(args.sdg_targets)

    # Load ES scores
    es_scores = defaultdict(lambda: defaultdict(float))  # entity → es_code → score
    es_names = {}
    with open(args.scores) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            entity = row['entity_id']
            es_code = row['es_code']
            es_scores[entity][es_code] = float(row['weighted_score'])
            es_names[es_code] = row.get('es_name', '')

    # Map ES scores → SDG targets
    # SDG score = sum(es_score × link_strength) for all linked ES codes
    sdg_target_scores = defaultdict(lambda: defaultdict(lambda: {
        'score': 0.0, 'contributing_es': [], 'goal_name': '', 'target_desc': '',
    }))

    for entity, scores in es_scores.items():
        for es_code, es_score in scores.items():
            if es_code not in crosswalk:
                continue
            for link in crosswalk[es_code]:
                target = link['sdg_target']
                contribution = es_score * link['link_strength']
                entry = sdg_target_scores[entity][target]
                entry['score'] += contribution
                entry['contributing_es'].append({
                    'es_code': es_code,
                    'es_name': es_names.get(es_code, ''),
                    'es_score': es_score,
                    'link_strength': link['link_strength'],
                    'mesh_es': link['mesh_es'],
                })
                entry['goal_name'] = link['goal_name']
                entry['target_desc'] = link.get('target_description', '')

    # Write per-target output
    target_cols = ['entity_id', 'sdg_goal', 'sdg_target', 'goal_name',
                   'target_description', 'sdg_score', 'contributing_es_codes',
                   'n_contributing_es']
    target_rows = []
    with open(args.output_targets, 'w') as f:
        f.write('\t'.join(target_cols) + '\n')
        for entity in sorted(sdg_target_scores):
            for target in sorted(sdg_target_scores[entity]):
                entry = sdg_target_scores[entity][target]
                goal = target.split('.')[0]
                es_list = ','.join(e['es_code'] for e in entry['contributing_es'])
                row_data = [
                    entity, goal, target, entry['goal_name'],
                    entry['target_desc'],
                    f"{entry['score']:.4f}", es_list,
                    str(len(entry['contributing_es'])),
                ]
                f.write('\t'.join(row_data) + '\n')
                target_rows.append({
                    'entity_id': entity,
                    'sdg_goal': goal,
                    'sdg_target': target,
                    'score': entry['score'],
                })

    # Aggregate to goal level
    goal_scores = defaultdict(lambda: defaultdict(lambda: {'score': 0.0, 'n_targets': 0}))
    for entity in sdg_target_scores:
        for target, entry in sdg_target_scores[entity].items():
            goal = target.split('.')[0]
            goal_scores[entity][goal]['score'] += entry['score']
            goal_scores[entity][goal]['n_targets'] += 1
            goal_scores[entity][goal]['goal_name'] = entry['goal_name']

    goal_cols = ['entity_id', 'sdg_goal', 'goal_name', 'sdg_score',
                 'n_targets', 'mean_target_score']
    with open(args.output_goals, 'w') as f:
        f.write('\t'.join(goal_cols) + '\n')
        for entity in sorted(goal_scores):
            for goal in sorted(goal_scores[entity], key=int):
                entry = goal_scores[entity][goal]
                mean = entry['score'] / entry['n_targets'] if entry['n_targets'] else 0
                f.write('\t'.join([
                    entity, goal, entry['goal_name'],
                    f"{entry['score']:.4f}",
                    str(entry['n_targets']),
                    f"{mean:.4f}",
                ]) + '\n')

    # Build viz JSON
    viz = {
        'targets': [],
        'goals': [],
        'entity_ids': sorted(sdg_target_scores.keys()),
    }

    for entity in sorted(sdg_target_scores):
        for target in sorted(sdg_target_scores[entity]):
            entry = sdg_target_scores[entity][target]
            viz['targets'].append({
                'entity_id': entity,
                'sdg_target': target,
                'sdg_goal': target.split('.')[0],
                'goal_name': entry['goal_name'],
                'score': round(entry['score'], 4),
                'contributing_es': [
                    {'es_code': e['es_code'], 'mesh_es': e['mesh_es'],
                     'contribution': round(e['es_score'] * e['link_strength'], 4)}
                    for e in entry['contributing_es']
                ],
            })

    for entity in sorted(goal_scores):
        for goal in sorted(goal_scores[entity], key=int):
            entry = goal_scores[entity][goal]
            viz['goals'].append({
                'entity_id': entity,
                'sdg_goal': int(goal),
                'goal_name': entry['goal_name'],
                'score': round(entry['score'], 4),
            })

    with open(args.output_json, 'w') as f:
        json.dump(viz, f, indent=2)

    # Stats
    total_entities = len(sdg_target_scores)
    total_goals = len(set(t.split('.')[0] for e in sdg_target_scores for t in sdg_target_scores[e]))
    total_targets = len(set(t for e in sdg_target_scores for t in sdg_target_scores[e]))
    print(f'SDG mapping: {total_entities} entities → {total_goals} SDG goals, '
          f'{total_targets} targets', file=sys.stderr)


if __name__ == '__main__':
    main()
