#!/usr/bin/env python3
"""
Compare two benchmark runs using bootstrap CI overlap.

Decision rule (from AGENT.md §4):
  - PROMOTE:      candidate CI lower > baseline CI upper  (real improvement)
  - INVESTIGATE:  CIs overlap but point estimates differ by > 0.02
  - DISCARD:      CIs overlap with point estimates within 0.02  (noise)

Also prints per-patient delta distribution to catch cases where
average improves but some patients tank.

Usage:
    python benchmark/compare_runs.py \
        --baseline experiments/baseline_2026-04-07/ \
        --candidate experiments/EXP-001/ \
        --metric recall_at_20
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_metric(run_dir, metric_name):
    """Load a metric JSON from a run directory."""
    path = Path(run_dir) / f'{metric_name}.json'
    if not path.exists():
        print(f"ERROR: {path} not found", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def load_per_patient(run_dir):
    """Load per-patient TSV from a run directory."""
    path = Path(run_dir) / 'per_patient.tsv'
    if not path.exists():
        print(f"ERROR: {path} not found", file=sys.stderr)
        sys.exit(1)
    return pd.read_csv(path, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="Compare two benchmark runs")
    parser.add_argument('--baseline', required=True, help='Path to baseline run directory')
    parser.add_argument('--candidate', required=True, help='Path to candidate run directory')
    parser.add_argument(
        '--metric', default='recall_at_20',
        choices=['recall_at_20', 'recall_at_50'],
        help='Metric to compare (default: recall_at_20)'
    )
    args = parser.parse_args()

    baseline = load_metric(args.baseline, args.metric)
    candidate = load_metric(args.candidate, args.metric)

    b_point = baseline['point_estimate']
    b_lo = baseline['ci_lower']
    b_hi = baseline['ci_upper']
    c_point = candidate['point_estimate']
    c_lo = candidate['ci_lower']
    c_hi = candidate['ci_upper']

    delta = c_point - b_point

    print(f"\n{'='*60}")
    print(f"  COMPARISON: {args.metric}")
    print(f"{'='*60}")
    print(f"\n  Baseline:  {b_point:.4f}  [{b_lo:.4f}, {b_hi:.4f}]")
    print(f"  Candidate: {c_point:.4f}  [{c_lo:.4f}, {c_hi:.4f}]")
    print(f"  Delta:     {delta:+.4f}")

    # Decision
    if c_lo > b_hi:
        decision = "PROMOTE"
        reason = "Candidate CI lower bound > baseline CI upper bound"
    elif abs(delta) > 0.02 and not (c_lo > b_hi):
        decision = "INVESTIGATE"
        reason = f"CIs overlap but point delta ({delta:+.4f}) > 0.02"
    else:
        decision = "DISCARD"
        reason = f"CIs overlap, point delta ({delta:+.4f}) within noise"

    print(f"\n  Decision:  ** {decision} **")
    print(f"  Reason:    {reason}")

    # Per-patient deltas
    metric_col = 'recall_at_20' if args.metric == 'recall_at_20' else 'recall_at_50'
    try:
        b_pp = load_per_patient(args.baseline)
        c_pp = load_per_patient(args.candidate)

        merged = b_pp.merge(
            c_pp, on='patient', suffixes=('_base', '_cand')
        )

        col_base = f'{metric_col}_base'
        col_cand = f'{metric_col}_cand'

        if col_base in merged.columns and col_cand in merged.columns:
            merged['delta'] = merged[col_cand] - merged[col_base]

            print(f"\n{'='*60}")
            print(f"  PER-PATIENT DELTAS ({len(merged)} patients)")
            print(f"{'='*60}")

            improved = (merged['delta'] > 0).sum()
            unchanged = (merged['delta'] == 0).sum()
            degraded = (merged['delta'] < 0).sum()
            tanked = (merged['delta'] < -0.5).sum()

            print(f"\n  Improved:  {improved}")
            print(f"  Unchanged: {unchanged}")
            print(f"  Degraded:  {degraded}")
            if tanked > 0:
                print(f"  ** TANKED (delta < -0.5): {tanked} patients **")

            print(f"\n  Delta distribution:")
            print(f"    min:    {merged['delta'].min():+.4f}")
            print(f"    25th:   {merged['delta'].quantile(0.25):+.4f}")
            print(f"    median: {merged['delta'].median():+.4f}")
            print(f"    75th:   {merged['delta'].quantile(0.75):+.4f}")
            print(f"    max:    {merged['delta'].max():+.4f}")

            # Show worst-hit patients
            worst = merged.nsmallest(5, 'delta')
            if not worst.empty and worst['delta'].iloc[0] < 0:
                print(f"\n  Worst-hit patients:")
                for _, row in worst.iterrows():
                    if row['delta'] < 0:
                        print(f"    Patient {row['patient']}: "
                              f"{row[col_base]:.3f} → {row[col_cand]:.3f} "
                              f"({row['delta']:+.4f})")

    except Exception as e:
        print(f"\n  (Could not compute per-patient deltas: {e})")

    print()
    return decision


if __name__ == '__main__':
    main()
