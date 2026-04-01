process MARK_DUPLICATES {
    tag "${patient_id} - ${sample_type}"
    label 'process_medium'
    container 'neoantigen/preprocessing:0.1.0'

    input:
    tuple val(patient_id), val(sample_type), path(bam), path(bai)

    output:
    tuple val(patient_id), val(sample_type), path("${patient_id}_${sample_type}.dedup.bam"), path("${patient_id}_${sample_type}.dedup.bam.bai"), emit: bam
    tuple val(patient_id), val(sample_type), path("dedup_metrics.txt"), emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${patient_id}_${sample_type}.dedup.bam \\
        -M dedup_metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT
    """
}

process BQSR {
    tag "${patient_id} - ${sample_type}"
    label 'process_medium'
    container 'neoantigen/preprocessing:0.1.0'

    input:
    tuple val(patient_id), val(sample_type), path(bam), path(bai)
    path reference
    path dbsnp
    path known_indels

    output:
    tuple val(patient_id), val(sample_type), path("${patient_id}_${sample_type}.recal.bam"), path("${patient_id}_${sample_type}.recal.bam.bai"), emit: bam

    script:
    def known_sites = ""
    if (dbsnp)        known_sites += " --known-sites ${dbsnp}"
    if (known_indels) known_sites += " --known-sites ${known_indels}"
    """
    gatk BaseRecalibrator \\
        -I ${bam} \\
        -R ${reference} \\
        ${known_sites} \\
        -O recal_table.txt

    gatk ApplyBQSR \\
        -I ${bam} \\
        -R ${reference} \\
        --bqsr-recal-file recal_table.txt \\
        -O ${patient_id}_${sample_type}.recal.bam

    samtools index ${patient_id}_${sample_type}.recal.bam
    """
}
