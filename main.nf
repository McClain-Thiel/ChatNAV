#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Neoantigen Vaccine Pipeline
 * ============================
 * Personalized mRNA neoantigen vaccine design from genomic sequencing data.
 *
 * Entry points are auto-detected from available inputs:
 *   - FASTQ  → full pipeline (modules 1-11)
 *   - BAM    → skip alignment (modules 2-11)
 *   - VCF    → skip preprocessing (modules 3/5-11)
 *   - Candidates + binding predictions → scoring only (modules 8-11)
 */

// ── Module includes ──

// Modules 1+2: Preprocessing
include { ALIGN }             from './modules/preprocessing/align'
include { MARKDUP_BQSR }      from './modules/preprocessing/markdup_bqsr'
include { MUTECT2 }           from './modules/preprocessing/mutect2'
include { CNVKIT }            from './modules/preprocessing/cnv'
include { STAR_FUSION }       from './modules/preprocessing/star_fusion'
include { VEP_ANNOTATE }      from './modules/preprocessing/vep'
include { PARABRICKS_FQ2BAM; PARABRICKS_MUTECT } from './modules/preprocessing/parabricks'

// Module 3: HLA typing
include { HLA_HD }            from './modules/hla_typing/hla_hd'
include { OPTITYPE }          from './modules/hla_typing/optitype'

// Module 4: Expression
include { SALMON_QUANT }      from './modules/expression/salmon'
include { GTEX_FALLBACK }     from './modules/expression/gtex_fallback'

// Module 5: Clonality
include { PYCLONE_VI }        from './modules/clonality/pyclone_vi'
include { SHARED_NEOANTIGEN } from './modules/clonality/shared_neoantigen'

// Module 6: Candidate generation
include { PVACSEQ }           from './modules/candidate_gen/pvacseq'
include { PVACFUSE }          from './modules/candidate_gen/pvacfuse'

// Module 7: MHC binding
include { MHC_BINDING_SUITE } from './subworkflows/mhc_binding_suite'

// Modules 8-11: Scoring, ranking, design
include { DOWNSTREAM }        from './subworkflows/downstream'

// Subworkflows
include { PREPROCESSING }     from './subworkflows/preprocessing'
include { PREPROCESSING_GPU } from './subworkflows/preprocessing_gpu'


// ── Entry point detection ──
def detect_entry_point(params) {
    if (params.binding_predictions && params.candidates_meta) {
        log.info "Entry point: pre-computed binding predictions → Modules 8-11"
        return 8
    }
    if (params.candidates_fasta && params.candidates_meta) {
        log.info "Entry point: candidate peptides → Modules 7-11"
        return 7
    }
    if (params.vcf) {
        log.info "Entry point: somatic VCF → Modules 5-11 (skipping preprocessing)"
        return 5
    }
    if (params.tumor_bam) {
        log.info "Entry point: aligned BAM → Modules 2-11"
        return 2
    }
    if (params.tumor_fastq_1) {
        log.info "Entry point: raw FASTQ → full pipeline (Modules 1-11)"
        return 1
    }
    error "No valid input files provided. Supply tumor_fastq_1, tumor_bam, vcf, or candidates_fasta."
}


workflow {

    def entry_point = detect_entry_point(params)

    // ════════════════════════════════════════════════════════════
    //  MODULES 1+2: Preprocessing + Variant Calling
    // ════════════════════════════════════════════════════════════

    if (entry_point <= 2) {
        tumor_fastq  = entry_point == 1
            ? Channel.of(tuple(params.patient_id, file(params.tumor_fastq_1), file(params.tumor_fastq_2)))
            : Channel.empty()
        normal_fastq = entry_point == 1
            ? Channel.of(tuple(params.patient_id, file(params.normal_fastq_1), file(params.normal_fastq_2)))
            : Channel.empty()
        tumor_bam_ch = entry_point == 2
            ? Channel.of(tuple(params.patient_id, file(params.tumor_bam)))
            : Channel.empty()
        normal_bam_ch = (entry_point == 2 && params.normal_bam)
            ? Channel.of(tuple(params.patient_id, file(params.normal_bam)))
            : Channel.empty()

        reference = file(params.reference_genome)

        if (params.parabricks) {
            PREPROCESSING_GPU(
                tumor_fastq,
                normal_fastq,
                reference
            )
            vcf_ch          = PREPROCESSING_GPU.out.vcf
            cnv_ch          = PREPROCESSING_GPU.out.cnv
            fusion_ch       = PREPROCESSING_GPU.out.fusions
            normal_bam_out  = PREPROCESSING_GPU.out.normal_bam
        } else {
            PREPROCESSING(
                tumor_fastq,
                normal_fastq,
                tumor_bam_ch,
                normal_bam_ch,
                reference
            )
            vcf_ch          = PREPROCESSING.out.vcf
            cnv_ch          = PREPROCESSING.out.cnv
            fusion_ch       = PREPROCESSING.out.fusions
            normal_bam_out  = PREPROCESSING.out.normal_bam
        }
    } else {
        // VCF entry point — read from params
        vcf_ch = Channel.of(tuple(params.patient_id, file(params.vcf)))
        cnv_ch = params.cnv_segments
            ? Channel.of(tuple(params.patient_id, file(params.cnv_segments)))
            : Channel.of(tuple(params.patient_id, file('NO_CNV')))
        fusion_ch = params.fusions
            ? Channel.of(tuple(params.patient_id, file(params.fusions)))
            : Channel.empty()
        normal_bam_out = params.normal_bam
            ? Channel.of(tuple(params.patient_id, file(params.normal_bam)))
            : Channel.empty()
    }


    // ════════════════════════════════════════════════════════════
    //  MODULE 3: HLA Typing
    // ════════════════════════════════════════════════════════════

    if (params.hla_alleles) {
        hla_ch = Channel.of(tuple(params.patient_id, file(params.hla_alleles)))
    } else if (entry_point <= 3) {
        HLA_HD(normal_bam_out)
        hla_ch = HLA_HD.out.alleles
    } else {
        error "HLA alleles must be provided when starting from VCF or later entry points."
    }


    // ════════════════════════════════════════════════════════════
    //  MODULE 4: Expression Quantification
    // ════════════════════════════════════════════════════════════

    if (params.expression) {
        expr_ch = Channel.of(tuple(params.patient_id, file(params.expression)))
    } else if (params.rnaseq_fastq_1 && entry_point <= 4) {
        rnaseq_ch = Channel.of(tuple(
            params.patient_id,
            file(params.rnaseq_fastq_1),
            file(params.rnaseq_fastq_2)
        ))
        SALMON_QUANT(rnaseq_ch, file(params.salmon_index))
        expr_ch = SALMON_QUANT.out.tpm
    } else if (params.tumor_type) {
        GTEX_FALLBACK(
            Channel.of(tuple(params.patient_id, params.tumor_type)),
            file(params.gtex_priors_table)
        )
        expr_ch = GTEX_FALLBACK.out.tpm
    } else {
        error "No expression data, RNA-seq FASTQ, or tumor_type for GTEx fallback."
    }


    // ════════════════════════════════════════════════════════════
    //  MODULES 5-7: Clonality → Candidates → MHC Binding
    // ════════════════════════════════════════════════════════════

    if (entry_point <= 7) {

        if (entry_point <= 5) {
            // Module 5: Clonality
            PYCLONE_VI(vcf_ch.join(cnv_ch))
            SHARED_NEOANTIGEN(
                PYCLONE_VI.out.ccf,
                file(params.shared_neoantigens_table)
            )
            ccf_ch = SHARED_NEOANTIGEN.out.ccf_annotated
        }

        if (entry_point <= 6) {
            // Module 6: Candidate generation
            ccf_vcf = (entry_point <= 5)
                ? vcf_ch.join(ccf_ch)
                : vcf_ch

            PVACSEQ(ccf_vcf.join(hla_ch).join(expr_ch))
            candidates_ch = PVACSEQ.out.candidates

            if (fusion_ch) {
                PVACFUSE(fusion_ch.join(hla_ch))
                candidates_ch = candidates_ch.mix(PVACFUSE.out.candidates)
            }
        } else {
            candidates_ch = Channel.of(tuple(
                params.patient_id,
                file(params.candidates_fasta),
                file(params.candidates_meta)
            ))
        }

        if (entry_point <= 7) {
            // Module 7: MHC binding prediction
            MHC_BINDING_SUITE(candidates_ch, hla_ch)
            binding_ch = MHC_BINDING_SUITE.out.binding_predictions
        }

    } else {
        // Entry point 8: pre-computed binding predictions
        binding_ch = Channel.of(tuple(
            params.patient_id,
            file(params.binding_predictions)
        ))
        candidates_ch = Channel.of(tuple(
            params.patient_id,
            file(params.candidates_fasta),
            file(params.candidates_meta)
        ))
    }


    // ════════════════════════════════════════════════════════════
    //  MODULES 8-11: Scoring → Ranking → Polyepitope → mRNA
    // ════════════════════════════════════════════════════════════

    DOWNSTREAM(
        binding_ch,
        candidates_ch,
        hla_ch,
        expr_ch,
        file(params.scoring_weights),
        file(params.utr_5prime),
        file(params.utr_3prime),
        file(params.signal_peptides)
    )
}
