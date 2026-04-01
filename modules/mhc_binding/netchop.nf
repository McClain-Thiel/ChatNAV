process NETCHOP {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/netmhc:0.1.0'

    input:
    tuple val(patient_id), path(candidates_fasta)

    output:
    tuple val(patient_id), path("netchop_output.txt"), emit: predictions

    script:
    """
    netchop ${candidates_fasta} > netchop_output.txt 2>/dev/null || true
    """
}
