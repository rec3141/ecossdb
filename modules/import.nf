// IMPORT_ANNOTATIONS — Normalize diverse annotation inputs

process IMPORT_ANNOTATIONS {
    label 'process_low'
    label 'ecossdb_core'

    publishDir "${params.outdir}/annotations", mode: 'copy', enabled: !params.store_dir
    storeDir   params.store_dir ? "${params.store_dir}/ecossdb/annotations" : null

    input:
    path annotation_files
    val  input_format

    output:
    path 'normalized_annotations.tsv', emit: tsv

    script:
    def format_flag = input_format == 'auto' ? '--format auto' : "--format ${input_format}"
    """
    normalize_annotations.py \
        --input ${annotation_files} \
        ${format_flag} \
        --output normalized_annotations.tsv
    """
}
