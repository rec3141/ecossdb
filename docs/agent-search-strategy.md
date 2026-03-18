# Agent Search Strategy: Literature-Based Mapping Expansion

## Current State

- **346 mappings** across 343 unique KOs
- **3 CICES classes covered**: 2.3.6.1 (climate), 2.3.5.1 (freshwater quality), 2.1.1.1 (waste decomposition)
- **Sources**: KEGG-Decoder (79), FOAM (196), cross-validated (71)
- **Missing**: 22 regulation/maintenance CICES classes have zero mappings

## Priority Targets (Aquatic Focus)

### Tier 1 — High-value unmapped classes with clear microbial links

| CICES | Name | Search strategy |
|-------|------|----------------|
| `2.3.5.2` | Chemical quality of salt water | Marine N/S/C cycling — mirror 2.3.5.1 with marine-specific genes |
| `2.3.3.2` | Controlling disease | Pathogen indicators, antimicrobial genes, ARGs |
| `2.1.1.2` | Filtering/sequestering pollutants | Heavy metal resistance, biosorption, bioaccumulation |
| `2.3.4.2` | Soil organic matter maintenance | C-degradation enzymes (CAZymes), humic compound processing |
| `2.3.6.2` | Local air quality | VOC degradation, DMS metabolism |

### Tier 2 — Expand existing classes with missing pathways

| Class | Missing pathways | Search strategy |
|-------|-----------------|----------------|
| `2.3.5.1` | P-cycling, Fe-cycling, Mn-cycling | Phosphatases, siderophores, Mn oxidation |
| `2.3.6.1` | N₂O flux, aerobic CH₄ oxidation specifics | nosZ clade II, pXMO variants |
| `2.1.1.1` | Pharmaceutical degradation, microplastics | Emerging contaminant degradation genes |

### Tier 3 — New ES classes with indirect microbial links

| Class | Approach |
|-------|----------|
| `2.2.1.1` | Erosion control → biofilm/EPS production genes |
| `2.3.4.1` | Soil formation → mineral weathering, rock phosphate solubilization |
| `2.3.3.1` | Pest control → insecticidal toxins, chitinases, biocontrol |

## Search Strategy per Target

Each target follows a 4-phase workflow:

### Phase 1: Anchor Paper Discovery

Find 2–3 review/meta-analysis papers that establish the gene–function–ES link.

**PubMed search templates:**

```
# Template A: Pathway-to-ES link
"{pathway}" AND ("ecosystem service" OR "biogeochemical") AND review[pt]

# Template B: Gene-to-function
"{gene_name}" AND ("{KO_id}" OR "{enzyme_name}") AND (metagenom* OR "functional annotation")

# Template C: Aquatic-specific
("{pathway}" OR "{function}") AND (freshwater OR marine OR groundwater OR "drinking water") AND microbi*
```

**Example for phosphorus cycling → 2.3.5.1:**
```
"phosphorus cycling" AND ("ecosystem service" OR "water quality") AND review[pt]
"polyphosphate" AND (ppk OR ppx) AND metagenom*
"phosphatase" AND (phoA OR phoD OR phoX) AND (freshwater OR marine)
```

### Phase 2: Gene Catalog Extraction

From anchor papers, extract specific gene–function claims:

1. **Read full text** (via PMC when available)
2. **Extract claims** of form: "Gene X (KO:Kxxxxx) performs function Y which contributes to ES Z"
3. **Score confidence** based on:
   - Direct experimental evidence in the paper → 0.85–0.95
   - Referenced from another study → 0.70–0.85
   - Inferred from pathway context → 0.50–0.70
   - Computational prediction only → 0.30–0.50

### Phase 3: Expand via Related Articles

For each anchor paper, use `find_related_articles` to discover:
- Papers citing the anchor (newer evidence)
- Papers with similar gene–function claims (cross-validation)
- Organism-specific studies (confidence boost for specific environments)

### Phase 4: Validate and Stage

1. Format as TSV rows matching `schema.json`
2. Run `validate_mapping.py` against ontology + existing mappings
3. Stage to `db/mappings/proposed/`

## Detailed Search Queries by Target

### Target: Phosphorus Cycling → 2.3.5.1 (Freshwater Quality)

```
Query 1: "phosphorus removal" AND (polyphosphate OR "phosphate accumulating") AND KO
Query 2: phoA[Title] OR phoD[Title] OR phoX[Title] AND phosphatase AND metagenom*
Query 3: "phosphorus cycling" AND (freshwater OR wastewater) AND "functional gene"
```

Expected genes: ppk1 (K00937), ppk2 (K22468), ppx (K01514), phoA (K01077), phoD (K01093), phoX (K07093), pstSCAB (K02040,K02037,K02038,K02036), pitA (K16322)

### Target: Iron/Manganese Cycling → 2.3.5.1

```
Query 1: ("iron oxidation" OR "iron reduction") AND (cyc2 OR mtrABC OR omcS) AND microb*
Query 2: ("manganese oxidation" OR "manganese reduction") AND metagenom* AND "functional gene"
Query 3: siderophore AND biosynthesis AND KO AND (freshwater OR marine)
```

Expected genes: mtrA (K21572), mtrB (K21573), mtrC (K21574), cymA (K15916), omcS, cyc2, mntABC

### Target: Pathogen/ARG Indicators → 2.3.3.2 (Disease Control)

```
Query 1: "antimicrobial resistance" AND "metagenomics" AND "indicator" AND "water quality"
Query 2: ("pathogen indicator" OR "fecal indicator") AND gene AND (freshwater OR "drinking water")
Query 3: "antibiotic resistance gene" AND ecosystem AND (ARG OR resistome) AND review[pt]
```

Expected genes: intI1 (integrase), sul1/sul2 (sulfonamide resistance), tetM/tetW, blaTEM, vanA, qnrS, mecA — mapped as inhibitors of disease control ES

### Target: Heavy Metal Resistance → 2.1.1.2 (Sequestering Pollutants)

```
Query 1: ("heavy metal" OR "metal resistance") AND (merA OR arsC OR chrA OR copA) AND metagenom*
Query 2: "bioremediation" AND ("metal sequestration" OR biosorption) AND "functional gene"
Query 3: ("mercury reduction" OR "arsenic resistance" OR "chromium reduction") AND KO
```

Expected genes: merA (K00520), merB (K00221), arsC (K00537), arsB (K03325), chrA (K07240), copA (K17686), czcABC, cadA, zntA

### Target: DMS/VOC Metabolism → 2.3.6.2 (Air Quality)

```
Query 1: "dimethylsulfoniopropionate" AND (dddP OR dddD OR dddL OR dmdA) AND marine
Query 2: "volatile organic compound" AND degradation AND microb* AND "functional gene"
Query 3: DMSP AND "climate" AND metagenom*
```

Expected genes: dddP (K16953), dddD, dddL, dmdA (K17486), isoprene degradation genes

### Target: CAZymes → 2.3.4.2 (Soil Organic Matter)

```
Query 1: "carbohydrate active enzyme" AND "soil organic matter" AND decomposition AND metagenom*
Query 2: (GH OR PL OR CE OR AA) AND lignocellulose AND "ecosystem service"
Query 3: "organic matter degradation" AND CAZy AND (soil OR sediment) AND "functional gene"
```

Expected: GH families (cellulases, xylanases, chitinases), AA families (laccases, peroxidases), PL families (pectinases)

## Agent Session Template

An agent executing this strategy should follow this workflow:

```
1. SELECT target CICES class from priority list
2. RUN 2-3 PubMed searches using templates above
3. READ top 3-5 review papers (full text via PMC when available)
4. EXTRACT gene-function-ES claims with evidence citations
5. FIND related articles for cross-validation
6. FORMAT as mapping TSV rows:
   - es_code: target CICES code
   - es_name: from ontology
   - gene_id: KEGG KO ID (preferred) or Pfam/COG
   - confidence: based on evidence rubric
   - evidence_source: literature:PMID:XXXXXXXX
   - functional_role: producer|consumer|transformer|inhibitor
   - pathway_context: specific pathway name
   - notes: one-line reasoning chain
7. VALIDATE with validate_mapping.py
8. STAGE to db/mappings/proposed/agent_<timestamp>.tsv
```

## Expected Yield

| Target | Est. new mappings | Confidence range |
|--------|-------------------|-----------------|
| P-cycling → 2.3.5.1 | 15–25 | 0.70–0.90 |
| Fe/Mn cycling → 2.3.5.1 | 10–20 | 0.60–0.85 |
| Marine quality → 2.3.5.2 | 80–120 | 0.50–0.80 |
| Pathogen/ARG → 2.3.3.2 | 20–40 | 0.65–0.85 |
| Heavy metals → 2.1.1.2 | 15–25 | 0.70–0.90 |
| DMS/VOC → 2.3.6.2 | 8–15 | 0.60–0.80 |
| CAZymes → 2.3.4.2 | 30–50 | 0.50–0.75 |
| **Total** | **~180–300** | |

This would roughly double the mapping table from 346 → 500–650 mappings and expand from 3 to ~8–10 CICES classes.

## Quality Controls

1. **No hallucinated KOs**: Every gene_id must be verified against KEGG (agent should check KO exists)
2. **PMID required**: Every literature-sourced mapping must cite a specific PMID
3. **Cross-validation boost**: If two independent papers support the same gene→ES link, confidence += 0.10
4. **Duplicate check**: validate_mapping.py catches exact duplicates against existing table
5. **Human review**: All agent-proposed mappings go to `proposed/` staging, never directly to master
