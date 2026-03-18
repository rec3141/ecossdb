# Gene-to-ES Mapping Guide

## Overview

The mapping table (`db/mappings/es_gene_mapping.tsv`) is the core database that links gene identifiers to ecosystem services. Each row maps one gene to one ES category.

## Columns

| Column | Required | Description |
|--------|----------|-------------|
| `es_code` | Yes | CICES code (e.g., `2.3.5.1`) |
| `es_name` | Yes | Human-readable ES name |
| `gene_id` | Yes | Gene identifier (e.g., `K00370`) |
| `gene_id_type` | Yes | One of: `KO`, `COG`, `Pfam`, `EC`, `CAZy`, `custom_hmm` |
| `confidence` | Yes | Score 0-1 (see rubric below) |
| `evidence_source` | Yes | Source of evidence (see formats below) |
| `functional_role` | Yes | One of: `producer`, `consumer`, `transformer`, `inhibitor` |
| `pathway_context` | No | Metabolic pathway context |
| `notes` | No | Free-text notes |

## Confidence Scoring Rubric

| Range | Level | Criteria | Example |
|-------|-------|----------|---------|
| 0.90-1.00 | Direct evidence | Gene product directly performs ES function, experimentally validated | mcrA → methanogenesis |
| 0.75-0.89 | Strong inference | Gene is part of a well-characterized pathway with clear ES link | narG → denitrification → water quality |
| 0.50-0.74 | Pathway-level | Gene is in a pathway associated with ES, but specific role less clear | General N-cycle enzyme |
| 0.25-0.49 | Indirect | Gene is indirectly related to ES through regulatory or auxiliary function | Transcription factor for N-cycle |
| 0.00-0.24 | Speculative | Theoretical link based on domain architecture or homology only | Novel protein with denitrification-like domain |

## Functional Roles

- **Producer**: The gene product directly provides or enhances the ecosystem service. Example: pmoA (methane oxidation) is a _producer_ of climate regulation because it removes CH4.
- **Consumer**: The gene product uses or depletes a resource related to the ES. Example: assimilatory nitrate reductase _consumes_ nitrate (reducing water purification load).
- **Transformer**: The gene product converts substrates in a way that contributes to the ES pathway. Example: nitrification enzymes _transform_ NH3→NO3.
- **Inhibitor**: The gene product actively works against the ES. Example: mcrA (methanogenesis) is an _inhibitor_ of climate regulation because it produces CH4.

## Evidence Source Formats

- `kegg_decoder` — Derived from KEGG-Decoder module definitions
- `foam` — Derived from FOAM ontology
- `kegg_decoder+foam` — Both sources agree
- `literature:PMID:12345678` — Literature citation by PubMed ID
- `literature:DOI:10.1234/example` — Literature citation by DOI
- `agent:<agent_id>` — Added by an AI agent

## How to Add Mappings

1. Create a TSV file with the same columns as the master table
2. Validate: `python3 bin/validate_mapping.py --ontology db/ontology/cices_v5.2.tsv --proposed your_file.tsv`
3. Place in `db/mappings/proposed/` for review
4. Merge: `python3 bin/merge_proposed_mappings.py --master db/mappings/es_gene_mapping.tsv --proposed your_file.tsv --output db/mappings/es_gene_mapping.tsv --log merge_log.tsv`

## Example Rows

```tsv
2.3.5.1	Controlling the chemical quality of freshwater	K00370	KO	0.85	kegg_decoder+literature	producer	Denitrification (NO3→N2)	narG nitrate reductase
2.3.6.1	Regulating our global climate	K00399	KO	0.90	kegg_decoder+literature	inhibitor	Methanogenesis (CO2→CH4)	mcrA; produces CH4
2.3.6.1	Regulating our global climate	K10944	KO	0.90	kegg_decoder	producer	Methane oxidation (CH4→CO2)	pmoA; consumes CH4
```
