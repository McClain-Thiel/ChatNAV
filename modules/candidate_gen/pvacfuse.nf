process PVACFUSE {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/pvactools:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_06", mode: 'copy'

    input:
    tuple val(patient_id), path(fusions_tsv), path(hla_alleles)

    output:
    tuple val(patient_id), path("fusion_candidates.fasta"), path("fusion_candidates_meta.tsv"), emit: candidates

    script:
    """
    HLA_LIST=\$(cat ${hla_alleles} | tr '\\n' ',' | sed 's/,\$//')

    pvacfuse run \\
        ${fusions_tsv} \\
        ${patient_id} \\
        "\${HLA_LIST}" \\
        NetMHCpan \\
        pvacfuse_output \\
        -e1 8,9,10,11 \\
        --n-threads ${task.cpus} \\
        --keep-tmp-files || true

    # Convert to pipeline format
    python3 << 'PYEOF'
import pandas as pd
import glob

results = []
for tsv in glob.glob("pvacfuse_output/*.filtered.tsv") + glob.glob("pvacfuse_output/*.all_epitopes.tsv"):
    try:
        df = pd.read_csv(tsv, sep='\\t')
        results.append(df)
    except Exception:
        pass

if not results:
    pd.DataFrame(columns=['peptide_id']).to_csv("fusion_candidates_meta.tsv", sep='\\t', index=False)
    with open("fusion_candidates.fasta", 'w') as f:
        pass
    exit(0)

combined = pd.concat(results, ignore_index=True)

meta = pd.DataFrame({
    'peptide_id': combined.get('MT Epitope Seq', ''),
    'peptide_sequence': combined.get('MT Epitope Seq', ''),
    'source_mutation': combined.get('Gene Name', 'fusion'),
    'mutation_type': 'fusion',
    'gene': combined.get('Gene Name', ''),
    'mutation': 'fusion',
    'ccf': 1.0,
    'tpm': combined.get('Gene Expression', 0),
    'wildtype_peptide': '',
    'is_frameshift': False,
    'is_shared': False,
    'hla_allele': combined.get('HLA Allele', ''),
    'mhc_class': 'I',
})

meta.to_csv("fusion_candidates_meta.tsv", sep='\\t', index=False)
with open("fusion_candidates.fasta", 'w') as f:
    for _, row in meta.drop_duplicates('peptide_id').iterrows():
        f.write(f">{row['peptide_id']}\\n{row['peptide_sequence']}\\n")
PYEOF
    """
}
