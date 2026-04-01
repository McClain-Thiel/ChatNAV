process SHARED_NEOANTIGEN {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/pyclone:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_05", mode: 'copy'

    input:
    tuple val(patient_id), path(ccf_tsv)
    path shared_table

    output:
    tuple val(patient_id), path("clonality_ccf_annotated.tsv"), emit: ccf_annotated

    script:
    """
    annotate_shared.py \\
        --ccf ${ccf_tsv} \\
        --shared-table ${shared_table} \\
        --output clonality_ccf_annotated.tsv
    """
}
