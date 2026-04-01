process ALIGN {
    tag "${patient_id} - ${sample_type}"
    label 'process_high'
    container 'neoantigen/preprocessing:0.1.0'

    input:
    tuple val(patient_id), val(sample_type), path(reads_1), path(reads_2), path(reference)

    output:
    tuple val(patient_id), val(sample_type), path("${patient_id}_${sample_type}.sorted.bam"), path("${patient_id}_${sample_type}.sorted.bam.bai"), emit: bam

    script:
    def rg = "@RG\\tID:${patient_id}_${sample_type}\\tSM:${patient_id}_${sample_type}\\tPL:ILLUMINA\\tLB:lib1"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R '${rg}' \\
        ${reference} \\
        ${reads_1} ${reads_2} \\
    | samtools sort -@ ${task.cpus} -o ${patient_id}_${sample_type}.sorted.bam

    samtools index ${patient_id}_${sample_type}.sorted.bam
    """
}
