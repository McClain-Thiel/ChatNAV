process STAR_FUSION {
    tag "${patient_id}"
    label 'process_high'
    container 'trinityctat/starfusion:1.12.0'
    publishDir "${params.outdir}/${patient_id}/module_02", mode: 'copy'

    input:
    tuple val(patient_id), path(reads_1), path(reads_2)
    path genome_lib

    output:
    tuple val(patient_id), path("fusions.tsv"), emit: fusions

    when:
    reads_1 && reads_2

    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${genome_lib} \\
        --left_fq ${reads_1} \\
        --right_fq ${reads_2} \\
        --CPU ${task.cpus} \\
        --output_dir star_fusion_output

    # Convert STAR-Fusion output to pipeline format
    python3 << 'PYEOF'
import pandas as pd
import os

sf_file = "star_fusion_output/star-fusion.fusion_predictions.tsv"
if not os.path.exists(sf_file):
    pd.DataFrame(columns=['gene1','gene2','junction_seq','tpm']).to_csv(
        "fusions.tsv", sep='\\t', index=False)
    exit(0)

sf = pd.read_csv(sf_file, sep='\\t')
result = pd.DataFrame({
    'gene1': sf['#FusionName'].str.split('--').str[0],
    'gene2': sf['#FusionName'].str.split('--').str[1],
    'junction_seq': sf.get('JunctionSequence', ''),
    'tpm': sf.get('FFPM', 0),
    'spanning_frags': sf.get('SpanningFragCount', 0),
    'junction_reads': sf.get('JunctionReadCount', 0),
})
result.to_csv("fusions.tsv", sep='\\t', index=False)
print(f"STAR-Fusion: {len(result)} gene fusions detected")
PYEOF
    """
}
