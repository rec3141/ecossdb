// SCORE_ES — Confidence-weighted ES scoring with pathway completeness

process SCORE_ES {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/scores", mode: 'copy', enabled: !params.store_dir

    input:
    path catalog
    path mapping_table
    val  role_weights

    output:
    path 'es_scores.tsv',     emit: scores
    path 'es_confidence.tsv', emit: confidence

    script:
    """
    score_es.py \
        --catalog ${catalog} \
        --mapping ${mapping_table} \
        --role-weights '${role_weights}' \
        --output es_scores.tsv \
        --confidence-out es_confidence.tsv
    """
}
