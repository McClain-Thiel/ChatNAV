#!/usr/bin/env python3
"""
Module 5 helper: Annotate variants with shared neoantigen status.

Matches variants by gene + protein change against the shared neoantigen
reference table (KRAS, TP53, IDH1, BRAF, etc.).

Inputs:
    clonality_ccf.tsv           — CCF estimates from PyClone-VI
    shared_neoantigens.tsv      — reference lookup table

Output:
    clonality_ccf_annotated.tsv — with is_shared_neoantigen column added
"""

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Annotate shared neoantigens')
    parser.add_argument('--ccf', required=True, help='clonality_ccf.tsv from PyClone-VI')
    parser.add_argument('--shared-table', required=True, help='shared_neoantigens.tsv reference')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    ccf = pd.read_csv(args.ccf, sep='\t')
    shared = pd.read_csv(args.shared_table, sep='\t')

    # Build lookup set: (gene, protein_change)
    shared_set = set()
    for _, row in shared.iterrows():
        shared_set.add((row['gene'].upper(), row['protein_change'].upper()))

    # Annotate
    ccf['is_shared_neoantigen'] = False
    for idx, row in ccf.iterrows():
        gene = str(row.get('gene', '')).upper()
        pchange = str(row.get('protein_change', '')).upper()
        if (gene, pchange) in shared_set:
            ccf.at[idx, 'is_shared_neoantigen'] = True

    shared_count = ccf['is_shared_neoantigen'].sum()
    print(f"Found {shared_count} shared neoantigens out of {len(ccf)} variants")

    ccf.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
