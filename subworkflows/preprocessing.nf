/*
 * Preprocessing subworkflow: Modules 1+2 (CPU path)
 * FASTQ → BWA-MEM2 → MarkDuplicates → BQSR → Mutect2 → VEP → annotated VCF
 * Also runs CNVkit and optionally STAR-Fusion
 */

include { ALIGN }              from '../modules/preprocessing/align'
include { MARK_DUPLICATES; BQSR } from '../modules/preprocessing/markdup_bqsr'
include { MUTECT2 }            from '../modules/preprocessing/mutect2'
include { CNVKIT }             from '../modules/preprocessing/cnv'
include { STAR_FUSION }        from '../modules/preprocessing/star_fusion'
include { VEP_ANNOTATE }       from '../modules/preprocessing/vep'

workflow PREPROCESSING {

    take:
    tumor_fastq_ch    // tuple(patient_id, reads_1, reads_2)
    normal_fastq_ch   // tuple(patient_id, reads_1, reads_2)
    tumor_bam_ch      // tuple(patient_id, bam) — for BAM entry point
    normal_bam_ch     // tuple(patient_id, bam) — for BAM entry point
    reference         // path to reference genome FASTA

    main:

    // ── Alignment (if starting from FASTQ) ──
    tumor_align_input = tumor_fastq_ch.map { pid, r1, r2 ->
        tuple(pid, 'tumor', r1, r2, reference)
    }
    normal_align_input = normal_fastq_ch.map { pid, r1, r2 ->
        tuple(pid, 'normal', r1, r2, reference)
    }

    ALIGN(tumor_align_input.mix(normal_align_input))

    // Split aligned BAMs by sample type
    tumor_aligned = ALIGN.out.bam
        .filter { it[1] == 'tumor' }
        .map { pid, stype, bam, bai -> tuple(pid, bam, bai) }
        .mix(tumor_bam_ch.map { pid, bam -> tuple(pid, bam, file("${bam}.bai")) })

    normal_aligned = ALIGN.out.bam
        .filter { it[1] == 'normal' }
        .map { pid, stype, bam, bai -> tuple(pid, bam, bai) }
        .mix(normal_bam_ch.map { pid, bam -> tuple(pid, bam, file("${bam}.bai")) })

    // ── MarkDuplicates ──
    tumor_dedup_input = tumor_aligned.map { pid, bam, bai ->
        tuple(pid, 'tumor', bam, bai)
    }
    normal_dedup_input = normal_aligned.map { pid, bam, bai ->
        tuple(pid, 'normal', bam, bai)
    }

    MARK_DUPLICATES(tumor_dedup_input.mix(normal_dedup_input))

    // ── BQSR ──
    dbsnp = params.dbsnp ? file(params.dbsnp) : file('NO_FILE')
    known_indels = params.known_indels ? file(params.known_indels) : file('NO_FILE')

    BQSR(MARK_DUPLICATES.out.bam, reference, dbsnp, known_indels)

    // Split recalibrated BAMs
    tumor_recal = BQSR.out.bam
        .filter { it[1] == 'tumor' }
        .map { pid, stype, bam, bai -> tuple(pid, bam, bai) }

    normal_recal = BQSR.out.bam
        .filter { it[1] == 'normal' }
        .map { pid, stype, bam, bai -> tuple(pid, bam, bai) }

    // ── Variant Calling ──
    gnomad = params.gnomad ? file(params.gnomad) : file('NO_FILE')
    MUTECT2(tumor_recal, normal_recal, reference, gnomad)

    // ── CNV Calling ──
    cnv_input = tumor_recal.join(normal_recal)
    CNVKIT(cnv_input, reference)

    // ── VEP Annotation ──
    vep_cache = file(params.vep_cache)
    VEP_ANNOTATE(MUTECT2.out.vcf, vep_cache, reference)

    // ── Fusion Detection (optional, if RNA-seq provided) ──
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
    normal_bam = normal_recal
}
