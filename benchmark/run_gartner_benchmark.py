#!/usr/bin/env python3
"""
Benchmark the pipeline against Gartner et al. 2021 (Nature Cancer) ground truth.

Uses Table 3 (confirmed CD8+ mmps with peptide sequences + HLA) and
Table 11 (all screened nmers with immunogenicity labels) to compute:
  - PPV@5, PPV@10, PPV@20 (positive predictive value at top K)
  - Recall@20 (fraction of true immunogenic peptides in top 20)
  - AUC for immunogenicity scoring

This benchmarks Modules 8a (BigMHC-IM) and 8b (structural scoring)
since we already have the peptides and HLA alleles from the Gartner data.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'bin'))


def load_gartner_data():
    """Load and merge Gartner ground truth data."""
    xlsx = PROJECT_ROOT / 'benchmark' / 'gartner' / 'supplementary_table.xlsx'

    # Table 3: Confirmed immunogenic mmps (peptide + HLA + features)
    positives = pd.read_excel(xlsx, sheet_name='Supplemetary Table 3', header=1)
    positives = positives.dropna(subset=['Mutant minimal Peptide', 'HLA'])

    # Table 1: Patient HLA alleles
    patients = pd.read_excel(xlsx, sheet_name='Supplementary Table 1', header=1)

    print(f"Loaded {len(positives)} confirmed immunogenic peptide-HLA pairs from {positives['ID'].nunique()} patients")
    return positives, patients


def score_with_bigmhc(peptides: list[str], hla_alleles: list[str]) -> list[float]:
    """Score peptides with BigMHC-IM."""
    from score_immunogenicity import _run_bigmhc_inprocess
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        for pep, hla in zip(peptides, hla_alleles):
            f.write(f"{pep},{hla}\n")
        input_path = f.name

    _run_bigmhc_inprocess(input_path)

    out_df = pd.read_csv(input_path + '.prd')
    Path(input_path).unlink(missing_ok=True)
    Path(input_path + '.prd').unlink(missing_ok=True)

    return out_df['BigMHC_IM'].tolist()


def score_with_structural(peptides: list[str], wt_peptides: list[str]) -> list[float]:
    """Score peptides with Tier 1 structural scoring."""
    from structural_scoring import score_tier1_position, find_mutation_positions

    scores = []
    for mut_pep, wt_pep in zip(peptides, wt_peptides):
        if wt_pep and isinstance(wt_pep, str) and len(wt_pep) == len(mut_pep):
            mut_pos = find_mutation_positions(mut_pep, wt_pep)
        else:
            mut_pos = list(range(len(mut_pep)))

        score, _ = score_tier1_position(mut_pep, mut_pos, len(mut_pep))
        scores.append(score if score is not None else 0.5)
    return scores


def compute_ppv_at_k(y_true: np.ndarray, y_score: np.ndarray, k: int) -> float:
    """Positive predictive value at top K."""
    top_k_idx = np.argsort(y_score)[::-1][:k]
    return np.mean(y_true[top_k_idx])


def compute_recall_at_k(y_true: np.ndarray, y_score: np.ndarray, k: int) -> float:
    """Recall at top K."""
    n_pos = np.sum(y_true)
    if n_pos == 0:
        return 0.0
    top_k_idx = np.argsort(y_score)[::-1][:k]
    return np.sum(y_true[top_k_idx]) / n_pos


def main():
    positives, patients = load_gartner_data()

    # We'll benchmark on the confirmed mmps — score them all, then check
    # if our scoring separates immunogenic from random peptides.

    # For a proper benchmark, we need BOTH positive and negative peptides.
    # The Gartner Table 3 only has positives. Let's create a matched negative
    # set by generating decoy peptides with the same length + HLA but different
    # (wildtype) sequences. The wildtype peptides are the natural negatives.

    peptides = []
    hla_alleles = []
    wt_peptides = []
    labels = []

    for _, row in positives.iterrows():
        mut_pep = str(row['Mutant minimal Peptide']).strip()
        wt_pep = str(row.get('Wild type minimal peptied', '')).strip()
        hla = str(row['HLA']).strip()

        if len(mut_pep) < 8 or len(mut_pep) > 14:
            continue
        if not hla.startswith('HLA-'):
            continue

        # Add mutant (positive)
        peptides.append(mut_pep)
        hla_alleles.append(hla)
        wt_peptides.append(wt_pep if wt_pep and wt_pep != '-' and wt_pep != 'nan' else '')
        labels.append(1)

        # Add wildtype as negative (if different from mutant)
        if wt_pep and wt_pep != '-' and wt_pep != 'nan' and wt_pep != mut_pep and len(wt_pep) == len(mut_pep):
            peptides.append(wt_pep)
            hla_alleles.append(hla)
            wt_peptides.append(wt_pep)
            labels.append(0)

    print(f"\nBenchmark set: {len(peptides)} peptides ({sum(labels)} positive, {len(labels) - sum(labels)} negative)")

    # Score with BigMHC-IM
    print("\nScoring with BigMHC-IM...")
    bigmhc_scores = score_with_bigmhc(peptides, hla_alleles)

    # Score with structural tier 1
    print("Scoring with structural tier 1...")
    structural_scores = score_with_structural(peptides, wt_peptides)

    # Composite score (simple weighted average for benchmark)
    composite = [0.7 * b + 0.3 * s for b, s in zip(bigmhc_scores, structural_scores)]

    labels_arr = np.array(labels)

    # Compute metrics
    print(f"\n{'='*60}")
    print(f"  BENCHMARK RESULTS — Gartner NCI Ground Truth")
    print(f"{'='*60}")
    print(f"\n  {len(peptides)} peptides: {sum(labels)} immunogenic + {len(labels) - sum(labels)} non-immunogenic")

    for name, scores in [('BigMHC-IM', bigmhc_scores),
                          ('Structural (Tier 1)', structural_scores),
                          ('Composite (0.7*IM + 0.3*Struct)', composite)]:
        scores_arr = np.array(scores)

        # AUC
        from sklearn.metrics import roc_auc_score, average_precision_score
        try:
            auc = roc_auc_score(labels_arr, scores_arr)
            ap = average_precision_score(labels_arr, scores_arr)
        except Exception:
            auc, ap = 0, 0

        print(f"\n  {name}:")
        print(f"    AUC:    {auc:.3f}")
        print(f"    AP:     {ap:.3f}")
        for k in [5, 10, 20, 50]:
            if k <= len(labels):
                ppv = compute_ppv_at_k(labels_arr, scores_arr, k)
                recall = compute_recall_at_k(labels_arr, scores_arr, k)
                print(f"    PPV@{k:<3d} {ppv:.3f}  |  Recall@{k:<3d} {recall:.3f}")

    # Per-patient metrics (macro average)
    print(f"\n  Per-patient breakdown (patients with >= 1 positive):")
    df = pd.DataFrame({
        'peptide': peptides,
        'hla': hla_alleles,
        'bigmhc': bigmhc_scores,
        'label': labels,
    })
    # Map peptide back to patient via positives table
    pep_to_patient = {}
    for _, row in positives.iterrows():
        pep = str(row['Mutant minimal Peptide']).strip()
        wt = str(row.get('Wild type minimal peptied', '')).strip()
        pid = row['ID']
        pep_to_patient[pep] = pid
        if wt and wt != '-' and wt != 'nan':
            pep_to_patient[wt] = pid

    df['patient'] = df['peptide'].map(pep_to_patient)

    # Print top 10 scored peptides
    print(f"\n  Top 10 scored peptides:")
    df_sorted = df.sort_values('bigmhc', ascending=False)
    for _, row in df_sorted.head(10).iterrows():
        marker = '+' if row['label'] == 1 else '-'
        print(f"    [{marker}] {row['peptide']:<15s} {row['hla']:<15s} BigMHC={row['bigmhc']:.4f}")


if __name__ == '__main__':
    main()
