process RANK_AND_SELECT {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/scoring:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_09", mode: 'copy'

    input:
    tuple val(patient_id), path(immunogenicity_scores), path(candidates_meta), path(binding_predictions)
    path scoring_weights

    output:
    tuple val(patient_id), path("selected_neoantigens.tsv"), emit: selected

    script:
    def structural_flag = params.structural_scoring ? '--structural-scoring' : ''
    """
    rank_and_select.py \\
        --immunogenicity-scores ${immunogenicity_scores} \\
        --candidates-meta ${candidates_meta} \\
        --binding-predictions ${binding_predictions} \\
        --scoring-weights ${scoring_weights} \\
        --weight-profile ${params.weight_profile} \\
        --top-n ${params.top_n} \\
        ${structural_flag} \\
        --output selected_neoantigens.tsv
    """
}
