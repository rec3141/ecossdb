// PARSE_ONTOLOGY — Convert CICES xlsx or hierarchical TSV to standardized format

process PARSE_ONTOLOGY {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/ontology", mode: 'copy', enabled: !params.store_dir
    storeDir   params.store_dir ? "${params.store_dir}/ecossdb/ontology" : null

    input:
    path ontology_input

    output:
    path 'es_ontology.tsv',    emit: tsv
    path 'es_hierarchy.json',  emit: json

    script:
    def is_xlsx = ontology_input.name.endsWith('.xlsx')
    if (is_xlsx)
        """
        parse_cices.py \
            --input ${ontology_input} \
            --output-tsv es_ontology.tsv \
            --output-json es_hierarchy.json
        """
    else
        """
        parse_ontology_tsv.py \
            --input ${ontology_input} \
            --output-tsv es_ontology.tsv \
            --output-json es_hierarchy.json
        """
}
