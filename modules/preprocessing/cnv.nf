process CNVKIT {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/cnvkit:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_02", mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference

    output:
    tuple val(patient_id), path("cnv_segments.tsv"), emit: cnv

    script:
    """
    cnvkit.py batch \\
        ${tumor_bam} \\
        --normal ${normal_bam} \\
        --fasta ${reference} \\
        --output-dir cnvkit_output \\
        -p ${task.cpus}

    # Convert CNVkit segments to pipeline format
    python3 << 'PYEOF'
import pandas as pd
import glob

seg_files = glob.glob("cnvkit_output/*.cns")
if not seg_files:
    # Create empty output
    pd.DataFrame(columns=['chr','start','end','log2_ratio','total_cn','minor_cn']).to_csv(
        "cnv_segments.tsv", sep='\\t', index=False)
    exit(0)

seg = pd.read_csv(seg_files[0], sep='\\t')
result = pd.DataFrame({
    'chr': seg['chromosome'],
    'start': seg['start'],
    'end': seg['end'],
    'log2_ratio': seg['log2'],
    'total_cn': (2 ** (seg['log2'] + 1)).round().astype(int).clip(0),
    'minor_cn': 0,  # CNVkit basic doesn't call allele-specific CN
})
result.to_csv("cnv_segments.tsv", sep='\\t', index=False)
print(f"CNVkit: {len(result)} segments")
PYEOF
    """
}
