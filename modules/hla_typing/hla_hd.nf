process HLA_HD_EXTRACT {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/hla-hd:0.1.0'

    input:
    tuple val(patient_id), path(normal_bam)

    output:
    tuple val(patient_id), path("hla_reads_1.fastq.gz"), path("hla_reads_2.fastq.gz"), emit: reads

    script:
    """
    # Extract reads from HLA region (chr6:28,477,797-33,448,354 in GRCh38)
    samtools view -b ${normal_bam} chr6:28477797-33448354 > hla_region.bam
    samtools sort -n hla_region.bam -o hla_sorted.bam
    samtools fastq -1 hla_reads_1.fastq.gz -2 hla_reads_2.fastq.gz \\
        -0 /dev/null -s /dev/null hla_sorted.bam
    """
}

process HLA_HD_RUN {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/hla-hd:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_03", mode: 'copy'

    input:
    tuple val(patient_id), path(reads_1), path(reads_2)

    output:
    tuple val(patient_id), path("hla_alleles.txt"), emit: alleles

    script:
    """
    # Run HLA-HD
    hlahd.sh \\
        -t ${task.cpus} \\
        -m 100 \\
        -f /opt/hlahd/freq_data \\
        ${reads_1} ${reads_2} \\
        /opt/hlahd/HLA_gene.split.txt \\
        /opt/hlahd/dictionary \\
        ${patient_id} \\
        hlahd_output

    # Parse HLA-HD results into simple allele list
    python3 << 'PYEOF'
import os
import glob

alleles = []
result_dir = glob.glob("hlahd_output/*/result")
if result_dir:
    for gene_file in glob.glob(f"{result_dir[0]}/*.txt"):
        with open(gene_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\\t')
                for part in parts[1:]:
                    if part.startswith('HLA-') or part.startswith('HLA_'):
                        allele = part.replace('_', '-')
                        if '*' in allele and ':' in allele:
                            # Keep 4-digit resolution: HLA-A*02:01
                            fields = allele.split(':')
                            allele = ':'.join(fields[:2])
                            alleles.append(allele)

# Deduplicate and sort
alleles = sorted(set(alleles))

with open("hla_alleles.txt", 'w') as f:
    for a in alleles:
        f.write(a + '\\n')

print(f"HLA-HD called {len(alleles)} alleles: {', '.join(alleles)}")
PYEOF
    """
}

workflow HLA_HD {
    take:
    normal_bam_ch   // tuple(patient_id, normal.bam)

    main:
    HLA_HD_EXTRACT(normal_bam_ch)
    HLA_HD_RUN(HLA_HD_EXTRACT.out.reads)

    emit:
    alleles = HLA_HD_RUN.out.alleles
}
