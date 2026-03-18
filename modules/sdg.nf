// MAP_ES_TO_SDG — Map ecosystem service scores to UN SDG targets

process MAP_ES_TO_SDG {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/sdg", mode: 'copy', enabled: !params.store_dir

    input:
    path scores
    path crosswalk
    path sdg_targets

    output:
    path 'es_sdg_targets.tsv', emit: targets
    path 'es_sdg_goals.tsv',   emit: goals
    path 'es_sdg.json',        emit: json

    script:
    """
    map_es_to_sdg.py \
        --scores ${scores} \
        --crosswalk ${crosswalk} \
        --sdg-targets ${sdg_targets} \
        --output-targets es_sdg_targets.tsv \
        --output-goals es_sdg_goals.tsv \
        --output-json es_sdg.json
    """
}
