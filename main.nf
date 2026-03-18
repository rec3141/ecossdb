#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ECOSSDB — Ecosystem Services Database Pipeline
 *
 * Maps metagenomic functional annotations to ecosystem services (ES)
 * using CICES 5.2 or custom ontologies.
 *
 * Usage:
 *   nextflow run main.nf --annotations <path> [--contig2bin <path>] [--proteins <path>]
 *   nextflow run main.nf -profile test
 */

// --- Help message ---
if (params.help) {
    log.info """
    ===========================================
     ECOSSDB v0.1 — Ecosystem Services Pipeline
    ===========================================

    Usage:
      nextflow run main.nf --annotations <path> [options]

    Required:
      --annotations     Path to annotation file(s) or directory

    Optional:
      --input_format    Input format: auto|danaseq|kofamscan|eggnog|gff|ko_list|dbcan (default: auto)
      --contig2bin      Contig-to-bin mapping TSV (enables per-MAG profiles)
      --proteins        Protein FASTA for custom HMM search
      --es_ontology_raw Raw ontology file (CICES xlsx or custom TSV) to parse
      --es_ontology     Pre-parsed ontology TSV (default: bundled CICES 5.2)
      --es_mapping      Gene-to-ES mapping table (default: bundled)
      --es_hmm_db       Custom HMM profiles (optional)
      --es_hmm_meta     HMM metadata TSV (required if es_hmm_db set)
      --sample_id       Sample identifier (default: 'sample')
      --outdir          Output directory (default: 'results')
      --role_weights    Scoring role weights (default: producer:1.0,transformer:0.7,consumer:-0.3,inhibitor:-0.5)
      --sdg             Enable SDG (Sustainable Development Goals) mapping (default: true)
    """.stripIndent()
    exit 0
}

// --- Validate inputs ---
if (!params.annotations) {
    error "ERROR: --annotations is required. Use --help for usage."
}

// --- Include modules ---
include { PARSE_ONTOLOGY        } from './modules/ontology'
include { IMPORT_ANNOTATIONS    } from './modules/import'
include { MAP_ANNOTATIONS_TO_ES } from './modules/mapping'
include { CUSTOM_HMM_SEARCH     } from './modules/mapping'
include { AGGREGATE_CONTIGS     } from './modules/aggregation'
include { AGGREGATE_BINS        } from './modules/aggregation'
include { SCORE_ES              } from './modules/scoring'
include { EXPORT_VIZ            } from './modules/viz'
include { MAP_ES_TO_SDG         } from './modules/sdg'

// --- Main workflow ---
workflow {

    // 1. Ontology
    if (params.es_ontology_raw) {
        ch_ontology_input = Channel.fromPath(params.es_ontology_raw)
        PARSE_ONTOLOGY(ch_ontology_input)
        ch_ontology_tsv  = PARSE_ONTOLOGY.out.tsv
        ch_hierarchy_json = PARSE_ONTOLOGY.out.json
    } else {
        ch_ontology_tsv   = Channel.fromPath(params.es_ontology)
        ch_hierarchy_json = Channel.fromPath("${projectDir}/db/ontology/es_hierarchy.json")
    }

    // 2. Import annotations
    ch_annotations = Channel.fromPath(params.annotations)
    IMPORT_ANNOTATIONS(ch_annotations.collect(), params.input_format)

    // 3. Map annotations to ES
    ch_mapping = Channel.fromPath(params.es_mapping)
    MAP_ANNOTATIONS_TO_ES(
        IMPORT_ANNOTATIONS.out.tsv,
        ch_mapping,
        ch_ontology_tsv
    )

    // 4. Optional custom HMM search
    if (params.es_hmm_db && params.es_hmm_meta) {
        ch_proteins = Channel.fromPath(params.proteins)
        ch_hmm_db   = Channel.fromPath(params.es_hmm_db)
        ch_hmm_meta = Channel.fromPath(params.es_hmm_meta)
        CUSTOM_HMM_SEARCH(ch_proteins, ch_hmm_db, ch_hmm_meta)
        ch_hmm_hits = CUSTOM_HMM_SEARCH.out.hits
    } else {
        ch_hmm_hits = Channel.fromPath("${projectDir}/assets/NO_HMM_HITS")
    }

    // 5. Aggregate to contigs
    AGGREGATE_CONTIGS(
        MAP_ANNOTATIONS_TO_ES.out.hits,
        ch_hmm_hits
    )

    // 6. Aggregate to bins (if contig2bin provided)
    if (params.contig2bin) {
        ch_contig2bin = Channel.fromPath(params.contig2bin)
        AGGREGATE_BINS(
            AGGREGATE_CONTIGS.out.catalog,
            ch_contig2bin
        )
        ch_mag_profiles = AGGREGATE_BINS.out.mag_profiles
    }

    // 7. Score ES
    SCORE_ES(
        AGGREGATE_CONTIGS.out.catalog,
        ch_mapping,
        params.role_weights
    )

    // 8. Export visualization JSON (if bins available)
    if (params.contig2bin) {
        EXPORT_VIZ(
            SCORE_ES.out.scores,
            AGGREGATE_CONTIGS.out.catalog,
            ch_mag_profiles,
            ch_hierarchy_json
        )
    }

    // 9. SDG mapping (optional, enabled by default)
    if (params.sdg) {
        ch_sdg_crosswalk = Channel.fromPath("${projectDir}/db/ontology/sdg/cices_to_sdg.tsv")
        ch_sdg_targets   = Channel.fromPath("${projectDir}/db/ontology/sdg/sdg_targets.tsv")
        MAP_ES_TO_SDG(
            SCORE_ES.out.scores,
            ch_sdg_crosswalk,
            ch_sdg_targets
        )
    }
}
