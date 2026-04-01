process MERGE_BINDING_RESULTS {
    tag "${patient_id}"
    label 'process_low'
    container 'neoantigen/mhc-binding:0.1.0'
    publishDir "${params.outdir}/${patient_id}/module_07", mode: 'copy'

    input:
    tuple val(patient_id), path(mhcflurry_tsv), path(mhcnuggets_tsv)

    output:
    tuple val(patient_id), path("binding_predictions.tsv"), emit: merged

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    dfs = []

    # MHCflurry Class I predictions
    try:
        df_i = pd.read_csv("${mhcflurry_tsv}", sep='\\t')
        if not df_i.empty:
            dfs.append(df_i)
            print(f"MHCflurry Class I: {len(df_i)} predictions")
    except Exception:
        pass

    # MHCnuggets Class II predictions
    try:
        df_ii = pd.read_csv("${mhcnuggets_tsv}", sep='\\t')
        if not df_ii.empty:
            dfs.append(df_ii)
            print(f"MHCnuggets Class II: {len(df_ii)} predictions")
    except Exception:
        pass

    if dfs:
        merged = pd.concat(dfs, ignore_index=True)
    else:
        merged = pd.DataFrame()

    # Ensure all expected columns exist
    for col in ['stability_rank', 'cleavage_score', 'tap_score',
                'presentation_score', 'processing_score']:
        if col not in merged.columns:
            merged[col] = None

    merged.to_csv("binding_predictions.tsv", sep='\\t', index=False)
    print(f"Merged binding predictions: {len(merged)} total")
    """
}
