// MAP_ANNOTATIONS_TO_ES — Core lookup: gene annotations → ES categories

process MAP_ANNOTATIONS_TO_ES {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/mapping", mode: 'copy', enabled: !params.store_dir

    input:
    path annotations
    path mapping_table
    path ontology

    output:
    path 'es_gene_hits.tsv',    emit: hits
    path 'es_mapping_stats.tsv', emit: stats

    script:
    """
    map_to_es.py \
        --annotations ${annotations} \
        --mapping ${mapping_table} \
        --ontology ${ontology} \
        --output es_gene_hits.tsv \
        --stats es_mapping_stats.tsv
    """
}


// CUSTOM_HMM_SEARCH — Optional custom HMM-based ES detection

process CUSTOM_HMM_SEARCH {
    label 'process_medium'
    label 'ecossdb_hmm'

    publishDir "${params.outdir}/mapping", mode: 'copy', enabled: !params.store_dir

    input:
    path proteins
    path hmm_db
    path hmm_metadata

    output:
    path 'hmm_es_hits.tsv', emit: hits

    script:
    """
    hmmsearch \
        --domtblout hmmsearch_results.domtblout \
        --noali \
        --cpu ${task.cpus} \
        ${hmm_db} ${proteins}

    hmm_to_es.py \
        --domtblout hmmsearch_results.domtblout \
        --metadata ${hmm_metadata} \
        --output hmm_es_hits.tsv
    """
}
