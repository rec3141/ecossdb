# ECOSSDB API for AI Agents

Machine-readable specification for extending the gene-to-ES mapping table.

## Schema

The mapping table schema is at `db/mappings/schema.json` (JSON Schema draft 2020-12).

Key constraints:
- `gene_id_type` enum: `KO`, `COG`, `Pfam`, `EC`, `CAZy`, `custom_hmm`
- `functional_role` enum: `producer`, `consumer`, `transformer`, `inhibitor`
- `confidence`: number in [0, 1]
- `es_code`: must match pattern `^[0-9]+(\.[0-9]+)*$` and exist in ontology

## Agent Workflow

### Step 1: Read Schema
```bash
cat db/mappings/schema.json
```

### Step 2: Read Ontology (valid ES codes)
```bash
# Get all valid ES codes and names
cut -f1,7 db/ontology/cices_v5.2.tsv
```

### Step 3: Check Existing Mappings
```bash
# Avoid duplicates
cat db/mappings/es_gene_mapping.tsv
```

### Step 4: Produce New Rows
Create a TSV file with header matching the schema columns:
```
es_code	es_name	gene_id	gene_id_type	confidence	evidence_source	functional_role	pathway_context	notes
```

Evidence sources should use structured references:
- `literature:PMID:12345678`
- `literature:DOI:10.1234/example`
- `agent:<your_agent_id>`

### Step 5: Validate
```bash
python3 bin/validate_mapping.py \
    --ontology db/ontology/cices_v5.2.tsv \
    --proposed new_rows.tsv \
    --existing db/mappings/es_gene_mapping.tsv \
    --json-output validation_results.json
```

The validation script returns:
```json
{
  "total": 10,
  "passed": 9,
  "failed": 1,
  "results": [
    {"line": 2, "gene_id": "K00370", "es_code": "2.3.5.1", "status": "pass", "errors": []},
    {"line": 3, "gene_id": "K99999", "es_code": "9.9.9.9", "status": "fail", "errors": ["ES code 9.9.9.9 not found in ontology"]}
  ]
}
```

### Step 6: Stage for Review
```bash
# Use timestamp in filename
cp new_rows.tsv db/mappings/proposed/agent_$(date +%Y%m%d_%H%M%S).tsv
```

### Step 7: Merge (after human approval)
```bash
python3 bin/merge_proposed_mappings.py \
    --master db/mappings/es_gene_mapping.tsv \
    --proposed db/mappings/proposed/agent_*.tsv \
    --output db/mappings/es_gene_mapping.tsv \
    --log merge_log.tsv
```

## Example: Complete Agent Session

An agent adding 3 nitrogen-cycle mappings:

```tsv
es_code	es_name	gene_id	gene_id_type	confidence	evidence_source	functional_role	pathway_context	notes
2.3.5.1	Controlling the chemical quality of freshwater	K00366	KO	0.80	literature:PMID:33456789	transformer	Assimilatory nitrite reduction	nirB ferredoxin-nitrite reductase
2.3.5.1	Controlling the chemical quality of freshwater	K17877	KO	0.75	literature:DOI:10.1038/s41586-020-12345	producer	Comammox nitrification	HAO-like hydroxylamine oxidase in comammox
2.1.1.1	Decomposing wastes or polluting substances	K14579	KO	0.70	agent:claude-2024	producer	Aromatic ring cleavage	Gentisate 1,2-dioxygenase
```

## Priority ES Codes for Aquatic Systems

These CICES codes are most relevant for microbial contributions in aquatic ecosystems:

| Code | Name | Key pathways |
|------|------|-------------|
| `2.3.5.1` | Controlling the chemical quality of freshwater | N-cycling, S-cycling, nutrient removal |
| `2.3.5.2` | Controlling the chemical quality of salt water | Marine biogeochemistry |
| `2.3.6.1` | Regulating our global climate | C-fixation, methanogenesis, methane oxidation |
| `2.1.1.1` | Decomposing wastes or polluting substances | Xenobiotic degradation, hydrocarbon breakdown |
| `2.1.1.2` | Filtering wastes or sequestering pollutants | Biosorption, bioaccumulation |
| `2.3.3.2` | Controlling disease | Pathogen indicators, antimicrobial genes |
