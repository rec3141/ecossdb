# ECOSSDB — Ecosystem Services Database Pipeline

Maps metagenomic functional annotations to ecosystem services using CICES 5.2.

## Quick Start

### Standalone

```bash
nextflow run main.nf \
    --annotations annotations.tsv \
    --contig2bin contig2bin.tsv \
    --outdir results
```

### With danaseq

```groovy
// In danaseq main.nf:
include { ECOSSDB } from '../ecossdb/adapters/danaseq'
ECOSSDB(MERGE_ANNOTATIONS.out.tsv, ch_contig2bin, ch_proteins)
```

### Test Profile

```bash
nextflow run main.nf -profile test
```

## Input Formats

ECOSSDB accepts annotations from multiple sources (auto-detected):

| Format | Description |
|--------|-------------|
| `danaseq` | danaseq merged_annotations.tsv |
| `kofamscan` | KofamScan detail output |
| `eggnog` | eggNOG-mapper annotations |
| `gff` | GFF3 with functional annotations |
| `ko_list` | Simple KO list (one per line) |
| `dbcan` | dbCAN overview output |

## Outputs

| File | Description |
|------|-------------|
| `es_gene_catalog.tsv` | Gene-level ES assignments |
| `es_per_contig.tsv` | Per-contig ES rollup |
| `es_per_mag.tsv` | Per-MAG ES profiles |
| `es_scores.tsv` | Confidence-weighted ES scores |
| `es_matrix.tsv` | Cross-sample ES matrix |
| `ecosystem_services.json` | Visualization data |

## Mapping Table

The gene-to-ES mapping table (`db/mappings/es_gene_mapping.tsv`) is the core database. It maps gene identifiers (KO, COG, Pfam, EC, CAZy) to CICES ecosystem service categories with confidence scores and functional roles.

Current coverage: **346 gene-ES mappings** across 3 CICES classes, bootstrapped from KEGG-Decoder modules and FOAM ontology.

### Extending the Mapping

See [docs/mapping-guide.md](docs/mapping-guide.md) for the human guide and [docs/api-for-agents.md](docs/api-for-agents.md) for the AI agent specification.

## Scoring

ES scores are computed as: `gene_presence × confidence × role_weight × pathway_completeness`

Role weights:
- **Producer** (+1.0): directly provides the ES
- **Transformer** (+0.7): converts substrates contributing to ES
- **Consumer** (-0.3): uses/depletes the ES
- **Inhibitor** (-0.5): actively works against the ES

See [docs/scoring.md](docs/scoring.md) for details.

## Custom Ontologies

ECOSSDB uses CICES 5.2 by default but accepts any hierarchical ES classification in TSV format. See [docs/ontology-guide.md](docs/ontology-guide.md).
