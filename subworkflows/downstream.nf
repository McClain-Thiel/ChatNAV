/*
 * Downstream subworkflow: Modules 8-11
 * Scoring → Ranking → Polyepitope Design → mRNA Design
 */

include { IMMUNOGENICITY_SCORING } from '../modules/scoring/immunogenicity'
include { STRUCTURAL_SCORING }     from '../modules/scoring/structural'
include { RANK_AND_SELECT }        from '../modules/ranking/rank_and_select'
include { POLYEPITOPE_DESIGN }     from '../modules/polyepitope/design'
include { MRNA_DESIGN }            from '../modules/mrna_design/mrna'

workflow DOWNSTREAM {

    take:
    binding_ch        // tuple(patient_id, binding_predictions.tsv)
    candidates_ch     // tuple(patient_id, candidates.fasta, candidates_meta.tsv)
    hla_ch            // tuple(patient_id, hla_alleles.txt)
    expr_ch           // tuple(patient_id, expression_tpm.tsv)
    scoring_weights   // path to scoring_weights.yaml
    utr_5prime        // path to 5' UTR FASTA
    utr_3prime        // path to 3' UTR FASTA
    signal_peptides   // path to signal_peptides.fasta

    main:

    // Prepare input for immunogenicity scoring:
    // Need (patient_id, binding_predictions, candidates_meta)
    // candidates_ch is (patient_id, fasta, meta) — extract meta
    candidates_meta_ch = candidates_ch.map { patient_id, fasta, meta ->
        tuple(patient_id, meta)
    }

    scoring_input = binding_ch.join(candidates_meta_ch)
    // Now: (patient_id, binding_predictions, candidates_meta)

    // Prepare human proteome path (optional)
    human_proteome = params.human_proteome
        ? file(params.human_proteome)
        : file('NO_FILE')

    // ── Module 8: Immunogenicity Scoring ──
    IMMUNOGENICITY_SCORING(scoring_input, human_proteome)
    immuno_scores = IMMUNOGENICITY_SCORING.out.scores

    // Optional structural scoring
    if (params.structural_scoring) {
        hla_chains = file(params.hla_heavy_chains)
        STRUCTURAL_SCORING(immuno_scores, hla_chains)
        final_scores = STRUCTURAL_SCORING.out.scores
    } else {
        final_scores = immuno_scores
    }

    // ── Module 9: Ranking and Selection ──
    // Need (patient_id, immuno_scores, candidates_meta, binding_predictions)
    ranking_input = final_scores
        .join(candidates_meta_ch)
        .join(binding_ch)
    // Now: (patient_id, immuno_scores, candidates_meta, binding_predictions)

    RANK_AND_SELECT(ranking_input, scoring_weights)

    // ── Module 10: Polyepitope Design ──
    POLYEPITOPE_DESIGN(RANK_AND_SELECT.out.selected, signal_peptides)

    // ── Module 11: mRNA Design ──
    MRNA_DESIGN(POLYEPITOPE_DESIGN.out.fasta, utr_5prime, utr_3prime)

    emit:
    selected    = RANK_AND_SELECT.out.selected
    polyepitope = POLYEPITOPE_DESIGN.out.fasta
    design      = POLYEPITOPE_DESIGN.out.design
    mrna        = MRNA_DESIGN.out.fasta
    spec        = MRNA_DESIGN.out.spec
}
