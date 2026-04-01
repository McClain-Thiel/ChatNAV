/*
 * MHC Binding Suite subworkflow: Module 7
 *
 * MHCflurry 2.0 (Class I) and MHCnuggets (Class II) run in parallel,
 * then results are merged into a unified binding_predictions.tsv.
 *
 * All tools are open-source and pip-installable — no academic licenses needed.
 *   - MHCflurry 2.0: Apache 2.0 license, includes antigen processing prediction
 *   - MHCnuggets: MIT license, supports Class I + II (we use it for Class II only)
 */

include { MHCFLURRY_PREDICT }    from '../modules/mhc_binding/mhcflurry'
include { MHCNUGGETS_PREDICT }   from '../modules/mhc_binding/mhcnuggets'
include { MERGE_BINDING_RESULTS } from '../modules/mhc_binding/merge_binding'

workflow MHC_BINDING_SUITE {

    take:
    candidates_ch   // tuple(patient_id, candidates.fasta, candidates_meta.tsv)
    hla_ch          // tuple(patient_id, hla_alleles.txt)

    main:

    // Extract FASTA from candidates channel
    fasta_ch = candidates_ch.map { patient_id, fasta, meta ->
        tuple(patient_id, fasta)
    }

    // Combine FASTA with HLA alleles
    fasta_hla = fasta_ch.join(hla_ch)
    // Now: (patient_id, candidates.fasta, hla_alleles.txt)

    // Run Class I and Class II binding prediction in parallel
    MHCFLURRY_PREDICT(fasta_hla)     // MHC-I: binding + antigen processing
    MHCNUGGETS_PREDICT(fasta_hla)    // MHC-II: binding affinity

    // Merge results
    merge_input = MHCFLURRY_PREDICT.out.predictions
        .join(MHCNUGGETS_PREDICT.out.predictions)

    MERGE_BINDING_RESULTS(merge_input)

    emit:
    binding_predictions = MERGE_BINDING_RESULTS.out.merged
}
