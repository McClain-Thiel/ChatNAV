#!/usr/bin/env python3
"""
Create a patient-stratified dev / muller-val split of the Muller dataset.

Splits by patient ID so the same patient never appears in both halves.
Stratifies by whether the patient has any immunogenic (CD8) mutations
to keep the positive-rate balanced across splits.

Usage:
    python benchmark/make_splits.py \
        --input benchmark/muller/Mutation_data_org.txt \
        --dev-frac 0.7 \
        --seed 42 \
        --output benchmark/muller/splits.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Create patient-stratified dev/val split")
    parser.add_argument("--input", required=True, help="Path to Mutation_data_org.txt")
    parser.add_argument("--dev-frac", type=float, default=0.7, help="Fraction of patients for dev set")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--output", required=True, help="Output path for splits.json")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')

    # Only consider screened mutations (CD8 or negative)
    screened = df[df['response_type'].isin(['CD8', 'negative'])]

    # Per-patient: does the patient have any immunogenic mutations?
    patient_has_pos = (
        screened.groupby('patient')['response_type']
        .apply(lambda x: (x == 'CD8').any())
        .reset_index()
        .rename(columns={'response_type': 'has_immunogenic'})
    )

    # Stratified split: sample dev-frac from each stratum
    rng = np.random.RandomState(args.seed)
    dev_patients = []
    val_patients = []

    for _, group in patient_has_pos.groupby('has_immunogenic'):
        pids = group['patient'].values.tolist()
        rng.shuffle(pids)
        n_dev = int(round(len(pids) * args.dev_frac))
        dev_patients.extend(pids[:n_dev])
        val_patients.extend(pids[n_dev:])

    # Sort for reproducibility in diffs
    dev_patients = sorted(dev_patients)
    val_patients = sorted(val_patients)

    splits = {
        "seed": args.seed,
        "dev_frac": args.dev_frac,
        "dev": dev_patients,
        "muller_val": val_patients,
        "n_dev": len(dev_patients),
        "n_val": len(val_patients),
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(splits, f, indent=2)

    # Summary
    dev_screened = screened[screened['patient'].isin(dev_patients)]
    val_screened = screened[screened['patient'].isin(val_patients)]
    dev_pos = (dev_screened['response_type'] == 'CD8').sum()
    val_pos = (val_screened['response_type'] == 'CD8').sum()

    print(f"Split {len(dev_patients) + len(val_patients)} patients "
          f"(seed={args.seed}, dev_frac={args.dev_frac}):")
    print(f"  dev:        {len(dev_patients)} patients, "
          f"{len(dev_screened)} screened mutations, {dev_pos} immunogenic")
    print(f"  muller-val: {len(val_patients)} patients, "
          f"{len(val_screened)} screened mutations, {val_pos} immunogenic")
    print(f"  Positive rate — dev: {dev_pos/len(dev_screened):.3f}, "
          f"val: {val_pos/len(val_screened):.3f}")
    print(f"Wrote {out_path}")


if __name__ == '__main__':
    main()
