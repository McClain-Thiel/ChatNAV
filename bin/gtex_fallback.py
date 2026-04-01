#!/usr/bin/env python3
"""
Module 4 fallback: Generate expression estimates from GTEx/TCGA priors
when RNA-seq data is not available.

All values are flagged as imputed=true.

Inputs:
    gtex_median_tpm_by_type.tsv     — median TPM by gene and tumor type
    tumor_type                       — e.g. 'SKCM'

Output:
    expression_tpm.tsv
"""

import argparse
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='GTEx expression fallback')
    parser.add_argument('--gtex-table', required=True)
    parser.add_argument('--tumor-type', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    gtex = pd.read_csv(args.gtex_table, sep='\t')
    tumor_type = args.tumor_type.upper()

    available = [c for c in gtex.columns if c not in ('gene_id', 'gene_name')]
    if tumor_type not in gtex.columns:
        raise ValueError(
            f"Tumor type '{tumor_type}' not found in GTEx table. "
            f"Available: {available}"
        )
        gtex['tpm'] = gtex[tumor_type]

    result = gtex[['gene_id', 'gene_name', 'tpm']].copy()
    result['imputed'] = True
    result.columns = ['gene_id', 'gene_name', 'tpm', 'expression_imputed']

    result.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(result)} imputed expression values for {tumor_type} → {args.output}")


if __name__ == '__main__':
    main()
