process MHCFLURRY_PREDICT {
    tag "${patient_id}"
    label 'process_medium'
    container 'neoantigen/mhc-binding:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_07", mode: 'copy', pattern: 'mhcflurry_*.tsv'

    input:
    tuple val(patient_id), path(candidates_fasta), path(hla_alleles)

    output:
    tuple val(patient_id), path("mhcflurry_predictions.tsv"), emit: predictions

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from mhcflurry import Class1PresentationPredictor

    # Load predictor (includes binding + antigen processing)
    predictor = Class1PresentationPredictor.load()

    # Read HLA alleles (Class I only — MHCflurry doesn't do Class II)
    with open("${hla_alleles}") as f:
        alleles = [line.strip() for line in f if line.strip()]
    class_i_alleles = [a for a in alleles if not any(
        x in a for x in ['DRB', 'DQB', 'DPB', 'DRA', 'DQA', 'DPA']
    )]

    if not class_i_alleles:
        raise ValueError("No Class I HLA alleles found")

    # Read candidate peptides from FASTA
    peptides = []
    with open("${candidates_fasta}") as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>') and line:
                peptides.append(line)

    if not peptides:
        pd.DataFrame().to_csv("mhcflurry_predictions.tsv", sep='\\t', index=False)
        exit(0)

    # MHCflurry 2.0 presentation prediction
    # Includes binding affinity + antigen processing (cleavage + TAP)
    results = predictor.predict(
        peptides=peptides,
        alleles=class_i_alleles,
        verbose=1,
    )

    # Rename to pipeline schema
    output = pd.DataFrame({
        'peptide_id': results['peptide'],
        'hla_allele': results['best_allele'],
        'mhc_class': 'I',
        'binding_rank': results['presentation_percentile'] / 100.0,
        'binding_affinity_nm': results.get('affinity', None),
        'presentation_score': results['presentation_score'],
        'processing_score': results['processing_score'],
        'stability_rank': None,  # MHCflurry integrates this into presentation score
        'cleavage_score': results.get('processing_score', None),
        'tap_score': None,  # integrated into processing_score
    })

    output.to_csv("mhcflurry_predictions.tsv", sep='\\t', index=False)
    print(f"MHCflurry: {len(output)} Class I predictions for {len(class_i_alleles)} alleles")
    """
}
