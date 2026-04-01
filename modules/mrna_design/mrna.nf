process MRNA_DESIGN {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/design:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_11", mode: 'copy'

    input:
    tuple val(patient_id), path(polyepitope_faa)
    path utr_5prime
    path utr_3prime

    output:
    tuple val(patient_id), path("mrna_sequence.fasta"),  emit: fasta
    tuple val(patient_id), path("synthesis_spec.json"),  emit: spec

    script:
    """
    design_mrna.py \\
        --polyepitope ${polyepitope_faa} \\
        --utr-5prime ${utr_5prime} \\
        --utr-3prime ${utr_3prime} \\
        --poly-a-length ${params.poly_a_length} \\
        --patient-id ${patient_id} \\
        --output-fasta mrna_sequence.fasta \\
        --output-spec synthesis_spec.json
    """
}
