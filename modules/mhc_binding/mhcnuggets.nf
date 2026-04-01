process MHCNUGGETS_PREDICT {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/mhc-binding:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_07", mode: 'copy', pattern: 'mhcnuggets_*.tsv'

    input:
    tuple val(patient_id), path(candidates_fasta), path(hla_alleles)

    output:
    tuple val(patient_id), path("mhcnuggets_predictions.tsv"), emit: predictions

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from mhcnuggets.src.predict import predict as mhcnuggets_predict
    import tempfile
    import os

    # Read HLA alleles — Class II only
    with open("${hla_alleles}") as f:
        alleles = [line.strip() for line in f if line.strip()]
    class_ii_alleles = [a for a in alleles if any(
        x in a for x in ['DRB', 'DQB', 'DPB']
    )]

    # Read peptides
    peptides = []
    with open("${candidates_fasta}") as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>') and line:
                peptides.append(line)

    if not class_ii_alleles or not peptides:
        # No Class II alleles or peptides — write empty output
        pd.DataFrame(columns=[
            'peptide_id','hla_allele','mhc_class','binding_rank'
        ]).to_csv("mhcnuggets_predictions.tsv", sep='\\t', index=False)
        print("No Class II alleles found — skipping MHCnuggets")
        exit(0)

    # Write peptides to temp file for MHCnuggets
    peptide_file = "peptides_for_mhcnuggets.txt"
    with open(peptide_file, 'w') as f:
        for p in peptides:
            f.write(p + '\\n')

    all_results = []
    for allele in class_ii_alleles:
        # MHCnuggets expects allele format: HLA-DRB10101 (no * or :)
        allele_fmt = allele.replace('*', '').replace(':', '')

        output_file = f"mhcnuggets_{allele_fmt}.csv"
        try:
            mhcnuggets_predict(
                class_='II',
                peptides_path=peptide_file,
                mhc=allele_fmt,
                output=output_file,
            )

            if os.path.exists(output_file):
                df = pd.read_csv(output_file)
                df['hla_allele'] = allele
                all_results.append(df)
                print(f"  {allele}: {len(df)} predictions")
        except Exception as e:
            print(f"  WARNING: MHCnuggets failed for {allele}: {e}")

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        # MHCnuggets outputs: peptide, ic50
        # Convert IC50 to percentile rank (approximate: <500nM ≈ <2%)
        combined['binding_rank'] = combined['ic50'].apply(
            lambda x: min(x / 25000.0, 1.0)  # rough IC50 → rank conversion
        )
        output = pd.DataFrame({
            'peptide_id': combined['peptide'],
            'hla_allele': combined['hla_allele'],
            'mhc_class': 'II',
            'binding_rank': combined['binding_rank'],
            'binding_affinity_nm': combined['ic50'],
            'presentation_score': None,
            'processing_score': None,
            'stability_rank': None,
            'cleavage_score': None,
            'tap_score': None,
        })
    else:
        output = pd.DataFrame(columns=[
            'peptide_id','hla_allele','mhc_class','binding_rank'
        ])

    output.to_csv("mhcnuggets_predictions.tsv", sep='\\t', index=False)
    print(f"MHCnuggets: {len(output)} Class II predictions")
    """
}
