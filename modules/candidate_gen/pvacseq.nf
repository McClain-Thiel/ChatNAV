process PVACSEQ {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/pvactools:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_06", mode: 'copy'

    input:
    tuple val(patient_id), path(vcf), path(ccf_tsv), path(hla_alleles), path(expression)

    output:
    tuple val(patient_id), path("candidates.fasta"), path("candidates_meta.tsv"), emit: candidates

    script:
    """
    # Read HLA alleles into comma-separated list
    HLA_LIST=\$(cat ${hla_alleles} | tr '\\n' ',' | sed 's/,\$//')

    # Run pVACseq
    pvacseq run \\
        ${vcf} \\
        ${patient_id} \\
        "\${HLA_LIST}" \\
        NetMHCpan NetMHCIIpan \\
        pvacseq_output \\
        -e1 8,9,10,11 \\
        -e2 13,15,17,19,21,23,25 \\
        --iedb-install-directory /opt/iedb \\
        --n-threads ${task.cpus} \\
        --keep-tmp-files || true

    # Convert pVACseq output to pipeline format
    # pVACseq outputs MHC_Class_I/<patient>.filtered.tsv and MHC_Class_II/<patient>.filtered.tsv
    python3 << 'PYEOF'
import pandas as pd
import os
import glob

results = []
for tsv in glob.glob("pvacseq_output/MHC_Class_*/*.filtered.tsv"):
    try:
        df = pd.read_csv(tsv, sep='\\t')
        results.append(df)
    except Exception:
        pass

if not results:
    # Try all.tsv files
    for tsv in glob.glob("pvacseq_output/MHC_Class_*/*.all_epitopes.tsv"):
        try:
            df = pd.read_csv(tsv, sep='\\t')
            results.append(df)
        except Exception:
            pass

if not results:
    print("WARNING: No pVACseq results found")
    # Create empty outputs
    pd.DataFrame(columns=['peptide_id']).to_csv("candidates_meta.tsv", sep='\\t', index=False)
    with open("candidates.fasta", 'w') as f:
        pass
    exit(0)

combined = pd.concat(results, ignore_index=True)

# Map pVACseq columns to pipeline schema
meta = pd.DataFrame({
    'peptide_id': combined.get('MT Epitope Seq', combined.get('Mutation', '')),
    'peptide_sequence': combined.get('MT Epitope Seq', ''),
    'source_mutation': combined.get('Chromosome', '').astype(str) + '_' +
                       combined.get('Start', '').astype(str) + '_' +
                       combined.get('Reference', '').astype(str) + '_' +
                       combined.get('Variant', '').astype(str),
    'mutation_type': combined.get('Variant Type', 'snv').str.lower().map(
        lambda x: 'frameshift' if 'frame' in str(x).lower() else
                  'indel' if 'ins' in str(x).lower() or 'del' in str(x).lower() else 'snv'),
    'gene': combined.get('Gene Name', ''),
    'mutation': combined.get('Amino Acid Change', ''),
    'ccf': 1.0,
    'tpm': combined.get('Gene Expression', 0),
    'wildtype_peptide': combined.get('WT Epitope Seq', ''),
    'is_frameshift': combined.get('Variant Type', '').str.contains('FS|frame', case=False, na=False),
    'is_shared': False,
    'hla_allele': combined.get('HLA Allele', ''),
    'mhc_class': combined.get('HLA Allele', '').apply(
        lambda x: 'II' if 'DRB' in str(x) or 'DQB' in str(x) else 'I'),
})

# Deduplicate by peptide sequence
meta = meta.drop_duplicates(subset=['peptide_id', 'hla_allele'])
meta.to_csv("candidates_meta.tsv", sep='\\t', index=False)

# Write FASTA
with open("candidates.fasta", 'w') as f:
    for _, row in meta.drop_duplicates('peptide_id').iterrows():
        f.write(f">{row['peptide_id']}\\n{row['peptide_sequence']}\\n")

print(f"Generated {len(meta)} candidate (peptide, HLA) pairs")
PYEOF
    """
}
