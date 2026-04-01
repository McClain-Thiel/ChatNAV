process STRUCTURAL_SCORING {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/scoring:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_08", mode: 'copy'

    input:
    tuple val(patient_id), path(immunogenicity_scores)
    path hla_heavy_chains

    output:
    tuple val(patient_id), path("immunogenicity_scores_structural.tsv"), emit: scores

    when:
    params.structural_scoring

    script:
    """
    structural_scoring.py \\
        --immunogenicity-scores ${immunogenicity_scores} \\
        --hla-heavy-chains ${hla_heavy_chains} \\
        --top-n 30 \\
        --output immunogenicity_scores_structural.tsv
    """
}
