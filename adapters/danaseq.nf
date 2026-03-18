#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * danaseq adapter — Maps danaseq pipeline outputs to ECOSSDB inputs
 *
 * This subworkflow takes danaseq output channels and feeds them
 * into the ECOSSDB pipeline. Include this in the danaseq main.nf
 * to add ecosystem service profiling.
 *
 * Usage in danaseq main.nf:
 *   include { ECOSSDB } from '../ecossdb/adapters/danaseq'
 *   ECOSSDB(MERGE_ANNOTATIONS.out.tsv, ch_contig2bin, ch_proteins)
 */

include { IMPORT_ANNOTATIONS    } from '../modules/import'
include { MAP_ANNOTATIONS_TO_ES } from '../modules/mapping'
include { AGGREGATE_CONTIGS     } from '../modules/aggregation'
include { AGGREGATE_BINS        } from '../modules/aggregation'
include { SCORE_ES              } from '../modules/scoring'
include { EXPORT_VIZ            } from '../modules/viz'

workflow ECOSSDB {
    take:
    ch_merged_annotations  // danaseq MERGE_ANNOTATIONS.out.tsv
    ch_contig2bin          // DAS_Tool or binning contig2bin
    ch_proteins            // Protein FASTA (optional, for HMM)

    main:
    // danaseq merged_annotations is already in normalized format
    IMPORT_ANNOTATIONS(ch_merged_annotations, 'danaseq')

    ch_ontology = Channel.fromPath(params.es_ontology ?: "${projectDir}/../ecossdb/db/ontology/cices_v5.2.tsv")
    ch_mapping  = Channel.fromPath(params.es_mapping  ?: "${projectDir}/../ecossdb/db/mappings/es_gene_mapping.tsv")
    ch_hierarchy = Channel.fromPath("${projectDir}/../ecossdb/db/ontology/es_hierarchy.json")

    MAP_ANNOTATIONS_TO_ES(
        IMPORT_ANNOTATIONS.out.tsv,
        ch_mapping,
        ch_ontology
    )

    // No HMM hits placeholder
    ch_hmm_placeholder = Channel.fromPath("${projectDir}/../ecossdb/assets/NO_HMM_HITS")

    AGGREGATE_CONTIGS(
        MAP_ANNOTATIONS_TO_ES.out.hits,
        ch_hmm_placeholder
    )

    AGGREGATE_BINS(
        AGGREGATE_CONTIGS.out.catalog,
        ch_contig2bin
    )

    SCORE_ES(
        AGGREGATE_CONTIGS.out.catalog,
        ch_mapping,
        params.role_weights ?: 'producer:1.0,transformer:0.7,consumer:-0.3,inhibitor:-0.5'
    )

    EXPORT_VIZ(
        SCORE_ES.out.scores,
        AGGREGATE_CONTIGS.out.catalog,
        AGGREGATE_BINS.out.mag_profiles,
        ch_hierarchy
    )

    emit:
    scores       = SCORE_ES.out.scores
    catalog      = AGGREGATE_CONTIGS.out.catalog
    mag_profiles = AGGREGATE_BINS.out.mag_profiles
    viz_json     = EXPORT_VIZ.out.json
}
