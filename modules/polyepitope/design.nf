process POLYEPITOPE_DESIGN {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/design:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_10", mode: 'copy'

    input:
    tuple val(patient_id), path(selected_neoantigens)
    path signal_peptides

    output:
    tuple val(patient_id), path("polyepitope.faa"),         emit: fasta
    tuple val(patient_id), path("polyepitope_design.tsv"),  emit: design

    script:
    """
    design_polyepitope.py \\
        --selected-neoantigens ${selected_neoantigens} \\
        --signal-peptides ${signal_peptides} \\
        --patient-id ${patient_id} \\
        --output-fasta polyepitope.faa \\
        --output-design polyepitope_design.tsv
    """
}
