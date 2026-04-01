process OPTITYPE {
    tag "${patient_id}"
    label 'process_medium'
    container 'fred2/optitype:1.3.5'
    publishDir "${params.outdir}/${patient_id}/module_03", mode: 'copy'

    input:
    tuple val(patient_id), path(tumor_bam)

    output:
    tuple val(patient_id), path("hla_alleles.txt"), emit: alleles

    script:
    """
    # Extract HLA region reads
    samtools view -b ${tumor_bam} chr6:28477797-33448354 > hla_region.bam
    samtools sort -n hla_region.bam -o hla_sorted.bam
    samtools fastq -1 hla_r1.fastq -2 hla_r2.fastq \\
        -0 /dev/null -s /dev/null hla_sorted.bam

    # OptiType (Class I only: HLA-A, B, C)
    OptiTypePipeline.py \\
        -i hla_r1.fastq hla_r2.fastq \\
        --dna \\
        -o optitype_output \\
        -p ${patient_id}

    # Parse OptiType TSV output
    python3 << 'PYEOF'
import pandas as pd
import glob

alleles = []
for tsv in glob.glob("optitype_output/*/*.tsv"):
    df = pd.read_csv(tsv, sep='\\t')
    for col in ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']:
        if col in df.columns:
            val = df[col].iloc[0]
            if pd.notna(val):
                allele = f"HLA-{col[0]}*{val}"
                alleles.append(allele)

alleles = sorted(set(alleles))
with open("hla_alleles.txt", 'w') as f:
    for a in alleles:
        f.write(a + '\\n')

print(f"OptiType called {len(alleles)} Class I alleles: {', '.join(alleles)}")
print("NOTE: OptiType does not call Class II alleles (DRB1, DQB1)")
PYEOF
    """
}
