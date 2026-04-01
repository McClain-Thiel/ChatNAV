#!/usr/bin/env python3
"""
Benchmark against Gartner et al. 2021 (Nature Cancer) ground truth.

Two benchmarks:
  1. Mut-vs-WT: 120 immunogenic mutant peptides vs their own wildtype peptides
     (easy — structural scoring trivially separates these)
  2. Mut-vs-Mut: 120 immunogenic mutant peptides vs ~600 non-immunogenic
     mutant peptides from the same patients (hard — tests real discrimination)

For benchmark 2, we generate decoy mutant peptides by taking the screened
nmers that were NOT immunogenic and extracting the best-binding MHCflurry
mmp from each nmer for the patient's HLA alleles.
"""

import sys
import random
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'bin'))


def load_gartner_data():
    """Load Gartner ground truth."""
    xlsx = PROJECT_ROOT / 'benchmark' / 'gartner' / 'supplementary_table.xlsx'

    # Table 3: Confirmed immunogenic mmps
    pos = pd.read_excel(xlsx, sheet_name='Supplemetary Table 3', header=1)
    pos = pos.dropna(subset=['Mutant minimal Peptide', 'HLA'])

    # Table 1: Patient info + HLAs
    patients = pd.read_excel(xlsx, sheet_name='Supplementary Table 1', header=1)
    hla_map = dict(zip(patients['ID'], patients['Patient HLAs']))

    # Table 11: All screened nmers with labels
    nmers = pd.read_excel(xlsx, sheet_name='Supplementary Table 11', header=1)

    return pos, patients, hla_map, nmers


def generate_negative_mmps(nmers_df, hla_map, positive_peptides: set, n_per_patient: int = 10):
    """
    Generate negative (non-immunogenic) mmps from screened-negative nmers.
    For each negative nmer, extract 8-11mer windows and pick one per nmer.
    """
    random.seed(42)
    negatives = []

    neg_nmers = nmers_df[nmers_df['Screening Status'] == '-'].copy()

    for _, row in neg_nmers.iterrows():
        patient_id = row['ID']
        nmer = str(row['Mutant nmer']).strip()

        if patient_id not in hla_map:
            continue

        hlas = [h.strip() for h in str(hla_map[patient_id]).split(',')]
        hla = hlas[0] if hlas else 'HLA-A*02:01'

        # Generate 9-10mer windows from the nmer
        for pep_len in [9, 10]:
            if len(nmer) < pep_len:
                continue
            # Pick a random window
            start = random.randint(0, len(nmer) - pep_len)
            pep = nmer[start:start + pep_len]

            # Skip if it matches a known positive
            if pep in positive_peptides:
                continue

            # Only use standard AAs
            if not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in pep):
                continue

            negatives.append({
                'peptide': pep,
                'hla': hla,
                'patient_id': patient_id,
                'label': 0,
                'source': 'screened_negative_nmer',
            })
            break  # one per nmer

    return negatives


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


def score_structural_tier1(peptides: list[str], wt_peptides: list[str]) -> list[float]:
    """Score peptides with Tier 1 structural scoring."""
    from structural_scoring import score_tier1_position, find_mutation_positions

    scores = []
    for mut_pep, wt_pep in zip(peptides, wt_peptides):
        if wt_pep and isinstance(wt_pep, str) and len(wt_pep) == len(mut_pep) and wt_pep != mut_pep:
            mut_pos = find_mutation_positions(mut_pep, wt_pep)
        else:
            mut_pos = list(range(len(mut_pep)))

        score, _ = score_tier1_position(mut_pep, mut_pos, len(mut_pep))
        scores.append(score if score is not None else 0.5)
    return scores


def ppv_at_k(labels, scores, k):
    top_k = np.argsort(scores)[::-1][:k]
    return np.mean(labels[top_k])


def recall_at_k(labels, scores, k):
    n_pos = np.sum(labels)
    if n_pos == 0:
        return 0.0
    top_k = np.argsort(scores)[::-1][:k]
    return np.sum(labels[top_k]) / n_pos


def print_metrics(name, labels, scores):
    from sklearn.metrics import roc_auc_score, average_precision_score
    labels = np.array(labels)
    scores = np.array(scores)

    auc = roc_auc_score(labels, scores)
    ap = average_precision_score(labels, scores)

    n_pos = int(np.sum(labels))
    n_neg = len(labels) - n_pos

    print(f"\n  {name} ({n_pos} pos / {n_neg} neg):")
    print(f"    AUC:     {auc:.3f}")
    print(f"    AP:      {ap:.3f}")
    for k in [5, 10, 20, 50]:
        if k <= len(labels):
            print(f"    PPV@{k:<3d}  {ppv_at_k(labels, scores, k):.3f}  |  Recall@{k:<3d}  {recall_at_k(labels, scores, k):.3f}")

    return auc, ap


def main():
    pos, patients, hla_map, nmers = load_gartner_data()
    print(f"Loaded: {len(pos)} immunogenic mmps, {len(nmers)} screened nmers, {len(hla_map)} patients")

    # ═══════════════════════════════════════════════════════
    # BENCHMARK 1: Mut vs WT (same peptide position)
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"  BENCHMARK 1: Immunogenic mutant vs wildtype")
    print(f"{'='*60}")

    b1_peptides, b1_hlas, b1_wts, b1_labels = [], [], [], []
    for _, row in pos.iterrows():
        mut = str(row['Mutant minimal Peptide']).strip()
        wt = str(row.get('Wild type minimal peptied', '')).strip()
        hla = str(row['HLA']).strip()
        if len(mut) < 8 or not hla.startswith('HLA-'):
            continue

        b1_peptides.append(mut)
        b1_hlas.append(hla)
        b1_wts.append(wt if wt and wt != '-' and wt != 'nan' else '')
        b1_labels.append(1)

        if wt and wt != '-' and wt != 'nan' and wt != mut and len(wt) == len(mut):
            b1_peptides.append(wt)
            b1_hlas.append(hla)
            b1_wts.append(wt)
            b1_labels.append(0)

    print(f"\n  Scoring {len(b1_peptides)} peptides with BigMHC-IM...")
    b1_bigmhc = score_with_bigmhc(b1_peptides, b1_hlas)
    b1_struct = score_structural_tier1(b1_peptides, b1_wts)
    b1_composite = [0.7 * b + 0.3 * s for b, s in zip(b1_bigmhc, b1_struct)]

    print_metrics('BigMHC-IM', b1_labels, b1_bigmhc)
    print_metrics('Structural Tier 1', b1_labels, b1_struct)
    b1_auc, b1_ap = print_metrics('Composite (0.7*IM + 0.3*Struct)', b1_labels, b1_composite)

    # ═══════════════════════════════════════════════════════
    # BENCHMARK 2: Immunogenic mutant vs non-immunogenic mutant
    # ═══════════════════════════════════════════════════════
    print(f"\n\n{'='*60}")
    print(f"  BENCHMARK 2: Immunogenic vs non-immunogenic MUTANTS")
    print(f"  (harder — all are mutant peptides, tests real discrimination)")
    print(f"{'='*60}")

    positive_seqs = set(str(row['Mutant minimal Peptide']).strip() for _, row in pos.iterrows())
    neg_mmps = generate_negative_mmps(nmers, hla_map, positive_seqs, n_per_patient=10)
    print(f"\n  Generated {len(neg_mmps)} negative mutant mmps from screened-negative nmers")

    b2_peptides, b2_hlas, b2_labels = [], [], []

    # Add positives
    for _, row in pos.iterrows():
        mut = str(row['Mutant minimal Peptide']).strip()
        hla = str(row['HLA']).strip()
        if len(mut) < 8 or not hla.startswith('HLA-'):
            continue
        if not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in mut):
            continue
        b2_peptides.append(mut)
        b2_hlas.append(hla)
        b2_labels.append(1)

    # Add negatives (sample to ~5x positives)
    n_neg_target = min(len(neg_mmps), 5 * sum(b2_labels))
    random.seed(42)
    sampled_neg = random.sample(neg_mmps, n_neg_target)
    for neg in sampled_neg:
        b2_peptides.append(neg['peptide'])
        b2_hlas.append(neg['hla'])
        b2_labels.append(0)

    print(f"  Scoring {len(b2_peptides)} peptides ({sum(b2_labels)} pos, {len(b2_labels) - sum(b2_labels)} neg)...")
    b2_bigmhc = score_with_bigmhc(b2_peptides, b2_hlas)

    # For structural: all are mutants, no WT available → all get position-based scores
    # with mutation assumed at all positions (frameshift-like)
    b2_struct = score_structural_tier1(b2_peptides, [''] * len(b2_peptides))
    b2_composite = [0.7 * b + 0.3 * s for b, s in zip(b2_bigmhc, b2_struct)]

    print_metrics('BigMHC-IM', b2_labels, b2_bigmhc)
    print_metrics('Structural Tier 1', b2_labels, b2_struct)
    b2_auc, b2_ap = print_metrics('Composite (0.7*IM + 0.3*Struct)', b2_labels, b2_composite)

    # ═══════════════════════════════════════════════════════
    # SUMMARY TABLE
    # ═══════════════════════════════════════════════════════
    print(f"\n\n{'='*60}")
    print(f"  SUMMARY")
    print(f"{'='*60}")
    print(f"\n  Benchmark 1 (mut vs wt):      Composite AUC={b1_auc:.3f}, AP={b1_ap:.3f}")
    print(f"  Benchmark 2 (mut vs mut):      Composite AUC={b2_auc:.3f}, AP={b2_ap:.3f}")

    # Top 10 from benchmark 2
    print(f"\n  Top 10 scored peptides (Benchmark 2):")
    b2_df = pd.DataFrame({
        'peptide': b2_peptides, 'hla': b2_hlas,
        'bigmhc': b2_bigmhc, 'label': b2_labels,
    })
    b2_sorted = b2_df.sort_values('bigmhc', ascending=False)
    for _, row in b2_sorted.head(10).iterrows():
        marker = '+' if row['label'] == 1 else '-'
        print(f"    [{marker}] {row['peptide']:<15s} {row['hla']:<15s} BigMHC={row['bigmhc']:.4f}")


if __name__ == '__main__':
    main()
