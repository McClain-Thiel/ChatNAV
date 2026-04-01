process VEP_ANNOTATE {
    tag "${patient_id}"
    label 'process_medium'
    container 'ensemblorg/ensembl-vep:release_110.1'
    publishDir "${params.outdir}/${patient_id}/module_02", mode: 'copy'

    input:
    tuple val(patient_id), path(vcf), path(tbi)
    path vep_cache
    path reference

    output:
    tuple val(patient_id), path("somatic_annotated.vcf.gz"), path("somatic_annotated.vcf.gz.tbi"), emit: vcf

    script:
    """
    vep \\
        --input_file ${vcf} \\
        --output_file somatic_annotated.vcf \\
        --format vcf \\
        --vcf \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --fasta ${reference} \\
        --assembly GRCh38 \\
        --offline \\
        --everything \\
        --pick \\
        --filter_common \\
        --fork ${task.cpus} \\
        --force_overwrite

    bgzip somatic_annotated.vcf
    tabix -p vcf somatic_annotated.vcf.gz
    """
}
