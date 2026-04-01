#!/usr/bin/env python3
"""
Module 7 helper: Merge outputs from NetMHCpan, NetMHCIIpan,
NetMHCstabpan, NetChop, and NetCTLpan into a unified TSV.

Each tool has an idiosyncratic output format. This script parses
each and joins on (peptide_id, hla_allele).

Output columns:
    peptide_id, hla_allele, mhc_class, binding_rank, binding_affinity_nm,
    stability_rank, cleavage_score, tap_score
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def parse_netmhcpan(path: str) -> pd.DataFrame:
    """Parse NetMHCpan 4.1 output (xls format)."""
    rows = []
    if not Path(path).exists():
        raise FileNotFoundError(f"Binding prediction file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            # Skip headers and separator lines
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            # NetMHCpan xls output is whitespace-delimited
            # Columns: Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM)
            fields = line.split()
            if len(fields) < 14:
                continue
            try:
                int(fields[0])  # Pos must be numeric
            except ValueError:
                continue

            rows.append({
                'peptide_id': fields[2],
                'hla_allele': fields[1].replace('*', '*'),
                'mhc_class': 'I',
                'binding_rank': float(fields[12]),  # %Rank_EL
                'binding_affinity_nm': float(fields[-1]) if len(fields) > 15 else None,
            })

    return pd.DataFrame(rows)


def parse_netmhciipan(path: str) -> pd.DataFrame:
    """Parse NetMHCIIpan 4.3 output."""
    rows = []
    if not Path(path).exists():
        raise FileNotFoundError(f"Binding prediction file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            fields = line.split()
            if len(fields) < 10:
                continue
            try:
                int(fields[0])
            except ValueError:
                continue

            rows.append({
                'peptide_id': fields[2],
                'hla_allele': fields[1],
                'mhc_class': 'II',
                'binding_rank': float(fields[8]),  # %Rank
                'binding_affinity_nm': float(fields[7]) if len(fields) > 7 else None,
            })

    return pd.DataFrame(rows)


def parse_netmhcstabpan(path: str) -> pd.DataFrame:
    """Parse NetMHCstabpan output — adds stability_rank."""
    rows = []
    if not Path(path).exists():
        raise FileNotFoundError(f"Binding prediction file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            fields = line.split()
            if len(fields) < 8:
                continue
            try:
                int(fields[0])
            except ValueError:
                continue

            rows.append({
                'peptide_id': fields[2],
                'hla_allele': fields[1],
                'stability_rank': float(fields[6]),  # %Rank
            })

    return pd.DataFrame(rows)


def parse_netchop(path: str) -> pd.DataFrame:
    """Parse NetChop output — proteasomal cleavage scores."""
    rows = []
    if not Path(path).exists():
        raise FileNotFoundError(f"Binding prediction file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            # NetChop: position aa score identity
            fields = line.split()
            if len(fields) < 4:
                continue
            try:
                float(fields[2])
            except ValueError:
                continue

            rows.append({
                'peptide_id': fields[3] if len(fields) > 3 else fields[0],
                'cleavage_score': float(fields[2]),
            })

    if not rows:
        return pd.DataFrame(rows)

    # Aggregate: take max cleavage score per peptide (C-terminal cleavage)
    df = pd.DataFrame(rows)
    return df.groupby('peptide_id').agg({'cleavage_score': 'max'}).reset_index()


def parse_netctlpan(path: str) -> pd.DataFrame:
    """Parse NetCTLpan output — TAP transport scores."""
    rows = []
    if not Path(path).exists():
        raise FileNotFoundError(f"Binding prediction file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            fields = line.split()
            if len(fields) < 8:
                continue
            try:
                int(fields[0])
            except ValueError:
                continue

            rows.append({
                'peptide_id': fields[2],
                'hla_allele': fields[1],
                'tap_score': float(fields[5]),
            })

    return pd.DataFrame(rows)


def parse_generic_tsv(path: str) -> pd.DataFrame:
    """Fallback: parse as a standard TSV (for pre-computed results)."""
    if not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t')


def main():
    parser = argparse.ArgumentParser(description='Merge MHC binding results')
    parser.add_argument('--netmhcpan', default=None, help='NetMHCpan 4.1 output')
    parser.add_argument('--netmhciipan', default=None, help='NetMHCIIpan 4.3 output')
    parser.add_argument('--netmhcstabpan', default=None, help='NetMHCstabpan output')
    parser.add_argument('--netchop', default=None, help='NetChop output')
    parser.add_argument('--netctlpan', default=None, help='NetCTLpan output')
    parser.add_argument('--precomputed', default=None, help='Pre-computed binding_predictions.tsv')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # If pre-computed results are provided, just pass through
    if args.precomputed:
        df = parse_generic_tsv(args.precomputed)
        df.to_csv(args.output, sep='\t', index=False)
        print(f"Using pre-computed binding predictions: {len(df)} rows")
        return

    # Parse each tool's output
    dfs = []
    if args.netmhcpan:
        mhc_i = parse_netmhcpan(args.netmhcpan)
        if not mhc_i.empty:
            dfs.append(mhc_i)
            print(f"NetMHCpan: {len(mhc_i)} predictions", file=sys.stderr)

    if args.netmhciipan:
        mhc_ii = parse_netmhciipan(args.netmhciipan)
        if not mhc_ii.empty:
            dfs.append(mhc_ii)
            print(f"NetMHCIIpan: {len(mhc_ii)} predictions", file=sys.stderr)

    if not dfs:
        print("ERROR: No binding predictions found", file=sys.stderr)
        pd.DataFrame().to_csv(args.output, sep='\t', index=False)
        return

    # Merge binding predictions
    merged = pd.concat(dfs, ignore_index=True)

    # Join stability scores
    if args.netmhcstabpan:
        stab = parse_netmhcstabpan(args.netmhcstabpan)
        if not stab.empty:
            merged = merged.merge(stab, on=['peptide_id', 'hla_allele'], how='left')
            print(f"NetMHCstabpan: {len(stab)} stability predictions", file=sys.stderr)

    # Join cleavage scores
    if args.netchop:
        chop = parse_netchop(args.netchop)
        if not chop.empty:
            merged = merged.merge(chop, on='peptide_id', how='left')
            print(f"NetChop: {len(chop)} cleavage scores", file=sys.stderr)

    # Join TAP scores
    if args.netctlpan:
        tap = parse_netctlpan(args.netctlpan)
        if not tap.empty:
            merged = merged.merge(tap, on=['peptide_id', 'hla_allele'], how='left')
            print(f"NetCTLpan: {len(tap)} TAP scores", file=sys.stderr)

    # Fill missing optional columns
    for col in ['stability_rank', 'cleavage_score', 'tap_score']:
        if col not in merged.columns:
            merged[col] = None

    merged.to_csv(args.output, sep='\t', index=False)
    print(f"Merged binding predictions: {len(merged)} rows → {args.output}")


if __name__ == '__main__':
    main()
