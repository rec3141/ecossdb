// EXPORT_VIZ — Generate visualization JSON for danaseq ES tab

process EXPORT_VIZ {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/viz", mode: 'copy', enabled: !params.store_dir

    input:
    path scores
    path catalog
    path mag_profiles
    path hierarchy_json

    output:
    path 'ecosystem_services.json',    emit: json
    path 'ecosystem_services.json.gz', emit: json_gz

    script:
    """
    es_to_json.py \
        --scores ${scores} \
        --catalog ${catalog} \
        --mag-profiles ${mag_profiles} \
        --hierarchy ${hierarchy_json} \
        --output ecosystem_services.json
    """
}
