process MUTECT2_CALL {
    tag "${patient_id} - ${interval}"
    label 'process_medium'
    container 'neoantigen/preprocessing:0.1.0'

    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(interval)
    path reference
    path gnomad

    output:
    tuple val(patient_id), path("${patient_id}_${interval}.vcf.gz"), path("${patient_id}_${interval}.vcf.gz.stats"), emit: vcf

    script:
    def gnomad_arg = gnomad ? "--germline-resource ${gnomad}" : ''
    """
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal \$(samtools view -H ${normal_bam} | grep '^@RG' | sed 's/.*SM:\\([^\\t]*\\).*/\\1/' | head -1) \\
        -L ${interval} \\
        ${gnomad_arg} \\
        -O ${patient_id}_${interval}.vcf.gz
    """
}

process MUTECT2_GATHER {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/preprocessing:0.1.0'

    input:
    tuple val(patient_id), path(vcfs), path(stats)

    output:
    tuple val(patient_id), path("${patient_id}.merged.vcf.gz"), path("${patient_id}.merged.stats"), emit: vcf

    script:
    def vcf_args = vcfs.collect { "-I ${it}" }.join(' ')
    def stats_args = stats.collect { "--stats ${it}" }.join(' ')
    """
    gatk GatherVcfs ${vcf_args} -O ${patient_id}.merged.vcf.gz

    gatk MergeMutectStats ${stats_args} -O ${patient_id}.merged.stats
    """
}

process MUTECT2_FILTER {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/preprocessing:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_02", mode: 'copy'

    input:
    tuple val(patient_id), path(vcf), path(stats)
    path reference

    output:
    tuple val(patient_id), path("${patient_id}.somatic.vcf.gz"), path("${patient_id}.somatic.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk FilterMutectCalls \\
        -R ${reference} \\
        -V ${vcf} \\
        --stats ${stats} \\
        -O ${patient_id}.filtered.vcf.gz

    # Keep only PASS variants with AF > 0.05 and depth > 10
    bcftools view \\
        -f PASS \\
        -i 'INFO/AF > 0.05 && INFO/DP > 10' \\
        ${patient_id}.filtered.vcf.gz \\
        -Oz -o ${patient_id}.somatic.vcf.gz

    tabix -p vcf ${patient_id}.somatic.vcf.gz
    """
}

workflow MUTECT2 {
    take:
    tumor_bam_ch    // tuple(patient_id, tumor.bam, tumor.bai)
    normal_bam_ch   // tuple(patient_id, normal.bam, normal.bai)
    reference
    gnomad

    main:
    // Scatter by chromosome
    chromosomes = Channel.of(
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
        'chr20','chr21','chr22','chrX','chrY'
    )

    // Combine BAMs with chromosomes for scatter
    scatter_input = tumor_bam_ch
        .join(normal_bam_ch)
        .combine(chromosomes)
    // Now: (patient_id, tumor.bam, tumor.bai, normal.bam, normal.bai, chr)

    MUTECT2_CALL(scatter_input, reference, gnomad)

    // Gather per patient
    gathered = MUTECT2_CALL.out.vcf
        .groupTuple(by: 0)
    // Now: (patient_id, [vcf1, vcf2, ...], [stats1, stats2, ...])

    MUTECT2_GATHER(gathered)
    MUTECT2_FILTER(MUTECT2_GATHER.out.vcf, reference)

    emit:
    vcf = MUTECT2_FILTER.out.vcf
}
