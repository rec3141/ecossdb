# Custom ES Ontology Guide

## Default: CICES 5.2

ECOSSDB ships with CICES V5.2 pre-parsed at `db/ontology/cices_v5.2.tsv`. To re-parse from the official xlsx:

```bash
python3 bin/parse_cices.py \
    --input /path/to/CICES_V5.2.xlsx \
    --output-tsv db/ontology/cices_v5.2.tsv \
    --output-json db/ontology/es_hierarchy.json
```

## Using a Custom Ontology

Any hierarchical classification in TSV format works. Provide it via `--es_ontology_raw`:

```bash
nextflow run main.nf \
    --annotations annotations.tsv \
    --es_ontology_raw my_ontology.tsv
```

### Required TSV Format

Minimum columns (additional columns are preserved but optional):

```
es_code	name	parent_code	level
1	Provisioning	-	section
1.1	Biomass	1	division
1.1.1	Wild plants	1.1	group
1.1.1.1	Edible wild plants	1.1.1	class
```

- `es_code`: Hierarchical code using dot notation
- `name`: Human-readable name
- `parent_code`: Code of the parent entry (`-` for root)
- `level`: Optional — inferred from code depth if omitted

### Full Format (matches CICES output)

```
es_code	level	section	division	group	class	name	description	parent_code
```

## Adding 5th-Digit CICES Extensions

CICES supports user-defined class types below the class level. Add subclass entries with 5-part codes:

```
es_code	name	parent_code	level
2.3.5.1.1	Nitrate removal	2.3.5.1	subclass
2.3.5.1.2	Phosphate removal	2.3.5.1	subclass
2.3.5.1.3	Heavy metal sequestration	2.3.5.1	subclass
```

These must also be added to the master ontology TSV and referenced in the mapping table.
