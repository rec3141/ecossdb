// AGGREGATE_CONTIGS — Merge hits, deduplicate, roll up to contigs

process AGGREGATE_CONTIGS {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/aggregation", mode: 'copy', enabled: !params.store_dir

    input:
    path gene_hits
    path hmm_hits

    output:
    path 'es_gene_catalog.tsv', emit: catalog
    path 'es_per_contig.tsv',   emit: contigs

    script:
    def hmm_flag = hmm_hits.name != 'NO_HMM_HITS' ? "--hmm-hits ${hmm_hits}" : ''
    """
    aggregate_contigs.py \
        --gene-hits ${gene_hits} \
        ${hmm_flag} \
        --output-catalog es_gene_catalog.tsv \
        --output-contigs es_per_contig.tsv
    """
}


// AGGREGATE_BINS — Roll up gene catalog to per-MAG profiles

process AGGREGATE_BINS {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/aggregation", mode: 'copy', enabled: !params.store_dir

    input:
    path catalog
    path contig2bin

    output:
    path 'es_per_mag.tsv', emit: mag_profiles

    script:
    """
    aggregate_bins.py \
        --catalog ${catalog} \
        --contig2bin ${contig2bin} \
        --output es_per_mag.tsv
    """
}


// AGGREGATE_SAMPLES — Cross-sample ES matrix

process AGGREGATE_SAMPLES {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/aggregation", mode: 'copy', enabled: !params.store_dir

    input:
    path mag_profiles
    val  sample_ids

    output:
    path 'es_matrix.tsv',  emit: matrix
    path 'es_summary.tsv', emit: summary

    script:
    def id_flag = sample_ids ? "--sample-ids ${sample_ids.join(' ')}" : ''
    """
    aggregate_samples.py \
        --mag-profiles ${mag_profiles} \
        ${id_flag} \
        --output-matrix es_matrix.tsv \
        --output-summary es_summary.tsv
    """
}
