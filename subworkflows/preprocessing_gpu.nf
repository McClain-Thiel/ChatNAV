/*
 * Preprocessing subworkflow: Modules 1+2 (GPU path / Parabricks)
 * FASTQ → Parabricks fq2bam → Parabricks mutectcaller → VEP
 */

include { PARABRICKS_FQ2BAM; PARABRICKS_MUTECT } from '../modules/preprocessing/parabricks'
include { CNVKIT }             from '../modules/preprocessing/cnv'
include { STAR_FUSION }        from '../modules/preprocessing/star_fusion'
include { VEP_ANNOTATE }       from '../modules/preprocessing/vep'

workflow PREPROCESSING_GPU {

    take:
    tumor_fastq_ch    // tuple(patient_id, reads_1, reads_2)
    normal_fastq_ch   // tuple(patient_id, reads_1, reads_2)
    reference         // path to reference genome FASTA

    main:

    // ── GPU-accelerated alignment + dedup + BQSR ──
    tumor_pb_input = tumor_fastq_ch.map { pid, r1, r2 ->
        tuple(pid, 'tumor', r1, r2, reference)
    }
    normal_pb_input = normal_fastq_ch.map { pid, r1, r2 ->
        tuple(pid, 'normal', r1, r2, reference)
    }

    PARABRICKS_FQ2BAM(tumor_pb_input.mix(normal_pb_input))

    tumor_bam = PARABRICKS_FQ2BAM.out.bam
        .filter { it[1] == 'tumor' }
        .map { pid, stype, bam, bai -> tuple(pid, bam, bai) }

    normal_bam = PARABRICKS_FQ2BAM.out.bam
        .filter { it[1] == 'normal' }
        .map { pid, stype, bam, bai -> tuple(pid, bam, bai) }

    // ── GPU-accelerated Mutect2 ──
    mutect_input = tumor_bam.join(normal_bam)
    PARABRICKS_MUTECT(mutect_input, reference)

    // ── VEP Annotation ──
    vep_cache = file(params.vep_cache)
    VEP_ANNOTATE(PARABRICKS_MUTECT.out.vcf, vep_cache, reference)

    // ── CNV (still CPU) ──
    cnv_input = tumor_bam.join(normal_bam)
    CNVKIT(cnv_input, reference)

    // ── Fusions ──
    fusions_ch = Channel.empty()
    if (params.rnaseq_fastq_1) {
        rnaseq_ch = Channel.of(tuple(
            params.patient_id,
            file(params.rnaseq_fastq_1),
            file(params.rnaseq_fastq_2)
        ))
        STAR_FUSION(rnaseq_ch, file(params.star_fusion_genome_lib ?: 'NO_FILE'))
        fusions_ch = STAR_FUSION.out.fusions
    }

    emit:
    vcf        = VEP_ANNOTATE.out.vcf
    cnv        = CNVKIT.out.cnv
    fusions    = fusions_ch
    normal_bam = normal_bam
}
