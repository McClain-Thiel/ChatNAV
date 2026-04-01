process GTEX_FALLBACK {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/scoring:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_04", mode: 'copy'

    input:
    tuple val(patient_id), val(tumor_type)
    path gtex_table

    output:
    tuple val(patient_id), path("expression_tpm.tsv"), emit: tpm

    script:
    """
    gtex_fallback.py \\
        --gtex-table ${gtex_table} \\
        --tumor-type ${tumor_type} \\
        --output expression_tpm.tsv
    """
}
