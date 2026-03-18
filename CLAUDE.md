# ECOSSDB Development Guide

## What is ECOSSDB?

ECOSSDB (Ecosystem Services Database) is a Nextflow pipeline that maps metagenomic functional annotations (KO, COG, Pfam, EC, CAZy) to ecosystem services using CICES 5.2 or custom ontologies. It takes functional annotation outputs from pipelines like danaseq, nf-core/mag, or standalone tools and produces ES profiles at gene, contig, MAG, and sample levels.

## Repository Structure

- `main.nf` — Pipeline entry point
- `modules/` — Nextflow process definitions (ontology, import, mapping, aggregation, scoring, viz)
- `bin/` — Python scripts (automatically on `$PATH` in Nextflow processes)
- `db/ontology/` — Parsed CICES 5.2 hierarchy (TSV + JSON)
- `db/mappings/` — Gene-to-ES mapping table, schema, proposed mappings staging
- `db/hmm/` — Custom HMM profiles and metadata
- `adapters/` — Integration subworkflows (danaseq)
- `test-data/` — Minimal test fixtures
- `tests/` — pytest suite
- `docs/` — Documentation

## Conventions

- Follows danaseq Nextflow patterns: `label` for resources, `conda` envs, `publishDir`/`storeDir` gating
- Python scripts use argparse, read/write TSV, print stats to stderr
- All intermediate data is plain TSV — git-diffable, spreadsheet-editable
- JSON only for viz output and ontology tree

## Process Labels

- `process_low`: 2 CPUs, 4 GB, 3h — most ECOSSDB processes
- `process_medium`: 8 CPUs, 16 GB, 12h — HMM search only
- `ecossdb_core`: conda env with Python, pandas, openpyxl, scipy
- `ecossdb_hmm`: conda env with HMMER3

## Adding a New Process

1. Create/edit the relevant module in `modules/`
2. Follow the pattern: label, conda, publishDir/storeDir, input/output/script
3. Add the Python script to `bin/` with `#!/usr/bin/env python3` shebang
4. Include the process in `main.nf`
5. Add test fixtures to `test-data/`

## Testing

```bash
# Run with test profile
nextflow run main.nf -profile test

# Run Python script tests
pytest tests/

# Test a single script
python3 bin/parse_cices.py --input /path/to/CICES.xlsx --output-tsv out.tsv --output-json out.json
```

## Mapping Table

The core data is `db/mappings/es_gene_mapping.tsv` — a plain TSV mapping gene IDs to CICES ES codes. See `db/mappings/schema.json` for the machine-readable schema and `docs/mapping-guide.md` for the human guide.

To extend the mapping table:
1. Add rows to a TSV file in `db/mappings/proposed/`
2. Validate: `python3 bin/validate_mapping.py --ontology db/ontology/cices_v5.2.tsv --proposed your_file.tsv --existing db/mappings/es_gene_mapping.tsv`
3. Merge: `python3 bin/merge_proposed_mappings.py --master db/mappings/es_gene_mapping.tsv --proposed db/mappings/proposed/*.tsv --output db/mappings/es_gene_mapping.tsv --log merge_log.tsv`

## Dependencies

All pipeline runtime scripts use **Python stdlib only** (csv, json, argparse, collections, re, gzip).
No pandas, numpy, scipy, or other third-party packages needed at runtime.

The only non-stdlib dependency is `openpyxl` in `parse_cices.py` — used once to parse the
CICES xlsx into TSV. The pre-parsed TSV is committed to the repo, so openpyxl is not needed
for normal pipeline execution.

### Standalone install
```bash
./install.sh          # Creates ecossdb-core conda env + verifies data files
./install.sh --skip-env  # Just verify data files (no conda)
```

### Within danaseq
No separate install needed. ECOSSDB is a git submodule initialized by danaseq's `install.sh`.
The ECOSSDB processes use danaseq's `dana-mag-classify` env (Python 3 + stdlib is sufficient).
Scripts are found via `beforeScript "export PATH=${ecossdb_bin}:$PATH"` in each process.

## Key Design Decisions

- CICES 5.2 as default but ontology-agnostic — any hierarchical TSV works
- Confidence-weighted scoring with functional role weights (producer, transformer, consumer, inhibitor)
- Gene-level granularity preserved through aggregation (contig → MAG → sample)
- AI agent-extensible: schema.json + validate_mapping.py + staging directory
- Pure stdlib Python — no heavy deps, runs in any Python 3.10+ env
