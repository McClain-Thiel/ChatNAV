process PARABRICKS_FQ2BAM {
    tag "${patient_id} - ${sample_type}"
    label 'process_gpu'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.1.0-1'

    input:
    tuple val(patient_id), val(sample_type), path(reads_1), path(reads_2), path(reference)

    output:
    tuple val(patient_id), val(sample_type), path("${patient_id}_${sample_type}.bam"), path("${patient_id}_${sample_type}.bam.bai"), emit: bam

    script:
    """
    pbrun fq2bam \\
        --ref ${reference} \\
        --in-fq ${reads_1} ${reads_2} \\
        --out-bam ${patient_id}_${sample_type}.bam \\
        --knownSites ${params.dbsnp} \\
        --num-gpus 1
    """
}

process PARABRICKS_MUTECT {
    tag "${patient_id}"
    label 'process_gpu'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.1.0-1'
    publishDir "${params.outdir}/${patient_id}/module_02", mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference

    output:
    tuple val(patient_id), path("${patient_id}.somatic.vcf.gz"), path("${patient_id}.somatic.vcf.gz.tbi"), emit: vcf

    script:
    """
    pbrun mutectcaller \\
        --ref ${reference} \\
        --in-tumor-bam ${tumor_bam} \\
        --in-normal-bam ${normal_bam} \\
        --out-vcf ${patient_id}.somatic.vcf.gz \\
        --num-gpus 1
    """
}
