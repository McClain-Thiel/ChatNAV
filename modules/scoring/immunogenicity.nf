process IMMUNOGENICITY_SCORING {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/scoring:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_08", mode: 'copy'

    input:
    tuple val(patient_id), path(binding_predictions), path(candidates_meta)
    path human_proteome

    output:
    tuple val(patient_id), path("immunogenicity_scores.tsv"), emit: scores

    script:
    def proteome_arg = human_proteome.name != 'NO_FILE' ? "--human-proteome ${human_proteome}" : ''
    """
    score_immunogenicity.py \\
        --binding-predictions ${binding_predictions} \\
        --candidates-meta ${candidates_meta} \\
        ${proteome_arg} \\
        --output immunogenicity_scores.tsv
    """
}
