process PREPARE_PYCLONE_INPUT {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/pyclone:0.1.0'

    input:
    tuple val(patient_id), path(vcf), path(cnv_segments)

    output:
    tuple val(patient_id), path("pyclone_input.tsv"), emit: input_tsv

    script:
    def cnv_arg = cnv_segments.name != 'NO_CNV' ? "--cnv-segments ${cnv_segments}" : ''
    """
    prepare_pyclone_input.py \\
        --vcf ${vcf} \\
        ${cnv_arg} \\
        --output pyclone_input.tsv
    """
}

process PYCLONE_VI_FIT {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/pyclone:0.1.0'

    input:
    tuple val(patient_id), path(pyclone_input)

    output:
    tuple val(patient_id), path("pyclone_results.tsv"), emit: results

    script:
    """
    pyclone-vi fit \\
        -i ${pyclone_input} \\
        -o pyclone_fit.h5 \\
        -c 40 \\
        -d beta-binomial \\
        -r 10

    pyclone-vi write-results-file \\
        -i pyclone_fit.h5 \\
        -o pyclone_results.tsv
    """
}

process FORMAT_CCF {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/pyclone:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_05", mode: 'copy'

    input:
    tuple val(patient_id), path(pyclone_results), path(pyclone_input)

    output:
    tuple val(patient_id), path("clonality_ccf.tsv"), emit: ccf

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    results = pd.read_csv("${pyclone_results}", sep='\\t')
    metadata = pd.read_csv("${pyclone_input}", sep='\\t')

    # PyClone-VI results: mutation_id, sample_id, cluster_id, cellular_prevalence, ...
    # Merge with original metadata for gene/consequence info
    merged = results.merge(
        metadata[['mutation_id', 'chr', 'position', 'gene', 'consequence', 'protein_change', 'vaf']],
        on='mutation_id', how='left'
    )

    # Rename to pipeline convention
    merged = merged.rename(columns={'cellular_prevalence': 'ccf'})
    merged['is_truncal'] = merged['ccf'] >= ${params.truncal_ccf}

    output_cols = ['mutation_id', 'chr', 'position', 'gene', 'consequence',
                   'protein_change', 'vaf', 'ccf', 'cluster_id', 'is_truncal']
    output_cols = [c for c in output_cols if c in merged.columns]
    merged[output_cols].to_csv("clonality_ccf.tsv", sep='\\t', index=False)
    print(f"CCF estimates: {len(merged)} variants, {merged['is_truncal'].sum()} truncal")
    """
}

workflow PYCLONE_VI {
    take:
    vcf_cnv_ch  // tuple(patient_id, vcf, cnv_segments)

    main:
    PREPARE_PYCLONE_INPUT(vcf_cnv_ch)
    PYCLONE_VI_FIT(PREPARE_PYCLONE_INPUT.out.input_tsv)
    FORMAT_CCF(
        PYCLONE_VI_FIT.out.results
            .join(PREPARE_PYCLONE_INPUT.out.input_tsv)
    )

    emit:
    ccf = FORMAT_CCF.out.ccf
}
