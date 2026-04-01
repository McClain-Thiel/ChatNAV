process SALMON_QUANT {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/salmon:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_04", mode: 'copy'

    input:
    tuple val(patient_id), path(reads_1), path(reads_2)
    path salmon_index

    output:
    tuple val(patient_id), path("expression_tpm.tsv"), emit: tpm

    script:
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${reads_1} \\
        -2 ${reads_2} \\
        -p ${task.cpus} \\
        --validateMappings \\
        -o salmon_output

    # Convert Salmon quant.sf to pipeline format
    python3 << 'PYEOF'
import pandas as pd

sf = pd.read_csv("salmon_output/quant.sf", sep='\\t')
result = pd.DataFrame({
    'gene_id': sf['Name'].str.split('|').str[0],
    'gene_name': sf['Name'].str.split('|').str[5] if '|' in sf['Name'].iloc[0] else sf['Name'],
    'tpm': sf['TPM'],
    'num_reads': sf['NumReads'],
    'expression_imputed': False,
})
result.to_csv("expression_tpm.tsv", sep='\\t', index=False)
print(f"Salmon quantification: {len(result)} transcripts, {result['tpm'].gt(1).sum()} with TPM > 1")
PYEOF
    """
}
