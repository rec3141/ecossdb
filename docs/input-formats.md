# Input Format Guide

## Auto-Detection

ECOSSDB auto-detects input format from file headers. Override with `--input_format <format>`.

## Supported Formats

### danaseq merged_annotations.tsv

The default format when running within danaseq. Tab-separated with header:

```
protein_id	contig_id	KO	COG	EC	Pfam	CAZy	description
contig001_1	contig001	K00370			PF00384		nitrate reductase alpha subunit
```

### KofamScan Detail Output

```
* contig001_1  K00370  nitrate reductase  1234.5  1.2e-100  ...
  contig001_2  K00371  nitrate reductase  456.7   3.4e-50   ...
```

Only `*`-marked (significant) hits are used.

### eggNOG-mapper Annotations

```
#query	seed_ortholog	evalue	score	...	KEGG_ko	EC	COG_category	...	PFAMs	...
contig001_1	...	...	...	...	ko:K00370	1.7.5.1	C	...	PF00384	...
```

### GFF3

```
##gff-version 3
contig001	prodigal	CDS	100	900	.	+	0	ID=contig001_1;KO=K00370;EC=1.7.5.1;product=nitrate reductase
```

### Simple KO List

One KO per line, or `protein_id<TAB>KO`:

```
K00370
K00371
K00399
```

or:

```
protein_001	K00370
protein_002	K00371
```

### dbCAN Overview

```
Gene ID	HMMER	DIAMOND	dbCAN_sub	#ofTools
contig001_1	GH13(1-400)	GH13	GH13_1	3
```

## Contig-to-Bin Mapping

Tab-separated, no header required:

```
contig001	bin_001
contig002	bin_002
```

Same format as DAS Tool `_DASTool_scaffolds2bin.txt` output.

## Standalone Usage with External Pipelines

### nf-core/mag

```bash
nextflow run main.nf \
    --annotations results/GenePrediction/prodigal/*/genes.gff \
    --input_format gff \
    --contig2bin results/GenomeBinning/DASTool/*/scaffolds2bin.txt
```

### MetaWRAP

```bash
nextflow run main.nf \
    --annotations metawrap_output/BIN_REFINEMENT/metawrap_bins.stats \
    --input_format ko_list \
    --contig2bin metawrap_output/BIN_REFINEMENT/contig2bin.tsv
```
