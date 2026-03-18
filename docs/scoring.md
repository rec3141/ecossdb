# ES Scoring Algorithm

## Overview

ECOSSDB computes ecosystem service scores using a composite formula:

```
ES_score = mean(gene_confidence × role_weight) × pathway_completeness
```

## Components

### Gene Confidence

Each gene-to-ES mapping has a confidence value (0-1) from the mapping table. This reflects the strength of evidence linking the gene to the ES.

### Functional Role Weights

| Role | Weight | Rationale |
|------|--------|-----------|
| Producer | +1.0 | Directly provides the ES |
| Transformer | +0.7 | Contributes to ES through metabolic transformation |
| Consumer | -0.3 | Depletes resources related to the ES |
| Inhibitor | -0.5 | Actively works against the ES |

Customizable via `--role_weights 'producer:1.0,transformer:0.7,consumer:-0.3,inhibitor:-0.5'`.

### Pathway Completeness

Analogous to KEGG module completeness. Genes are grouped by `pathway_context` within each ES code. Each unique pathway context is a "step"; completeness = fraction of steps with at least one gene present.

```
completeness = steps_present / total_steps
```

If no pathway information is available, completeness defaults to 1.0.

## Worked Example

Consider a MAG with the following genes mapped to ES code 2.3.5.1 (Water quality):

| Gene | Confidence | Role | Weight | Pathway |
|------|-----------|------|--------|---------|
| K00370 (narG) | 0.85 | producer | +1.0 | Denitrification |
| K00371 (narH) | 0.85 | producer | +1.0 | Denitrification |
| K00368 (nirK) | 0.85 | producer | +1.0 | Denitrification |
| K00376 (nosZ) | 0.85 | producer | +1.0 | Denitrification |
| K10944 (amoA) | 0.80 | transformer | +0.7 | Nitrification |

**Step 1**: Weighted contributions per gene
- narG: 0.85 × 1.0 = 0.85
- narH: 0.85 × 1.0 = 0.85
- nirK: 0.85 × 1.0 = 0.85
- nosZ: 0.85 × 1.0 = 0.85
- amoA: 0.80 × 0.7 = 0.56

**Step 2**: Mean weighted score
- mean = (0.85 + 0.85 + 0.85 + 0.85 + 0.56) / 5 = 0.792

**Step 3**: Pathway completeness
- Denitrification pathways present: 4 steps represented
- Nitrification pathways present: 1 step represented
- If denitrification has 4 total steps and nitrification has 3: completeness = 5/7 = 0.714

**Step 4**: Final score
- score = 0.792 × 0.714 = 0.565

## Output Files

- `es_scores.tsv`: Entity-level ES scores with completeness
- `es_confidence.tsv`: Per-gene detail showing confidence and role weights applied
