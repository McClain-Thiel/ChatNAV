#!/usr/bin/env python3
"""
Benchmark against Muller/Gfeller harmonized dataset (Immunity 2023).

131 patients, 48,306 mutations, 13,483 screened (213 immunogenic).
All signals available: expression (TPM), CCF, clonality, wildtype sequences,
pre-computed binding ranks (NetMHCpan, PRIME, MixMHCpred).

Generates ALL 8-11mer peptide windows per mutation, scores them all with
BigMHC-IM, and takes the best score per mutation — matching how the real
pipeline works.
"""

import sys
import math
import tempfile
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'bin'))

AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


def load_muller():
    df = pd.read_csv(PROJECT_ROOT / 'benchmark' / 'muller' / 'Mutation_data_org.txt', sep='\t')
    screened = df[df['response_type'].isin(['CD8', 'negative'])].copy()
    screened['label'] = (screened['response_type'] == 'CD8').astype(int)
    return screened


def generate_all_windows(mut_seq, wt_seq):
    """Generate all 8-11mer windows that contain the mutation."""
    windows = []
    if not mut_seq or mut_seq == 'nan' or len(mut_seq) < 8:
        return windows
    if not wt_seq or wt_seq == 'nan':
        wt_seq = ''

    for pep_len in [8, 9, 10, 11]:
        for start in range(len(mut_seq) - pep_len + 1):
            pep = mut_seq[start:start + pep_len]
            if not all(c in AMINO_ACIDS for c in pep):
                continue
            wt_pep = wt_seq[start:start + pep_len] if len(wt_seq) >= start + pep_len else ''
            # Only keep windows that actually contain the mutation
            if wt_pep and len(wt_pep) == len(pep) and wt_pep == pep:
                continue
            windows.append((pep, wt_pep))
    return windows


def run_bigmhc_batch(peptides, hla_alleles):
    from score_immunogenicity import _run_bigmhc_inprocess
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        for pep, hla in zip(peptides, hla_alleles):
            f.write(f"{pep},{hla}\n")
        path = f.name
    _run_bigmhc_inprocess(path)
    out = pd.read_csv(path + '.prd')
    Path(path).unlink(missing_ok=True)
    Path(path + '.prd').unlink(missing_ok=True)
    return dict(zip(out['pep'].tolist(), out['BigMHC_IM'].tolist()))


def compute_foreignness_batch(peptides, proteome_kmers):
    from score_immunogenicity import compute_foreignness_kmer
    results = {}
    for pep in peptides:
        if pep in results:
            continue
        try:
            results[pep] = compute_foreignness_kmer(pep, proteome_kmers)
        except ValueError:
            results[pep] = 0.5
    return results


def compute_structural_batch(peptides, wt_peptides):
    from structural_scoring import score_tier1_position, find_mutation_positions
    results = {}
    for mut, wt in zip(peptides, wt_peptides):
        if mut in results:
            continue
        if wt and isinstance(wt, str) and len(wt) == len(mut) and wt != mut:
            pos = find_mutation_positions(mut, wt)
        else:
            pos = list(range(len(mut)))
        s, _ = score_tier1_position(mut, pos, len(mut))
        results[mut] = s if s is not None else 0.5
    return results


def normalize_expression(tpm):
    if pd.isna(tpm) or tpm < 0:
        return 0.0
    return min(math.log10(tpm + 1) / 5.0, 1.0)


def normalize_agretopicity(agreto):
    if pd.isna(agreto):
        return 0.5
    return 1.0 / (1.0 + math.exp(-0.5 * agreto))


def ppv_at_k(labels, scores, k):
    idx = np.argsort(scores)[::-1][:k]
    return np.mean(labels[idx])


def recall_at_k(labels, scores, k):
    n_pos = np.sum(labels)
    if n_pos == 0:
        return 0.0
    idx = np.argsort(scores)[::-1][:k]
    return np.sum(labels[idx]) / n_pos


def auc_score(labels, scores):
    from sklearn.metrics import roc_auc_score
    try:
        return roc_auc_score(labels, scores)
    except ValueError:
        return 0.0


def print_metrics(name, labels, scores):
    labels = np.array(labels)
    scores = np.array(scores)
    a = auc_score(labels, scores)
    print(f"\n  {name}:")
    print(f"    AUC:     {a:.3f}")
    for k in [5, 10, 20, 50, 100]:
        if k <= len(labels):
            print(f"    PPV@{k:<4d} {ppv_at_k(labels, scores, k):.3f}  |  Recall@{k:<4d} {recall_at_k(labels, scores, k):.3f}")
    return a


def main():
    print("Loading Muller/Gfeller harmonized dataset...")
    df = load_muller()
    print(f"  {len(df)} screened mutations, {df['label'].sum()} immunogenic, {df['patient'].nunique()} patients")

    # ── Generate ALL peptide windows per mutation ──
    print("\n  Generating all 8-11mer peptide windows per mutation...")
    mutation_windows = {}  # mutation_idx → [(pep, wt_pep), ...]
    all_unique_peps = set()

    for idx, row in df.iterrows():
        windows = generate_all_windows(
            str(row.get('mutant_seq', '')),
            str(row.get('wt_seq', ''))
        )
        mutation_windows[idx] = windows
        for pep, _ in windows:
            all_unique_peps.add(pep)

    n_windows = sum(len(w) for w in mutation_windows.values())
    print(f"  {n_windows} total windows, {len(all_unique_peps)} unique peptides")

    # ── Score ALL unique peptides with BigMHC-IM ──
    unique_pep_list = list(all_unique_peps)
    hla_list = ['HLA-A*02:01'] * len(unique_pep_list)

    print(f"\n  Running BigMHC-IM on {len(unique_pep_list)} unique peptides...")
    bigmhc_lookup = run_bigmhc_batch(unique_pep_list, hla_list)
    print(f"  BigMHC scored {len(bigmhc_lookup)} peptides")

    # ── Foreignness for all unique peptides ──
    from score_immunogenicity import build_proteome_kmers
    proteome_path = PROJECT_ROOT / 'reference' / 'proteome' / 'human_proteome.fasta'
    print("  Building proteome k-mer index...")
    proteome_kmers = build_proteome_kmers(str(proteome_path))
    print(f"  Computing foreignness...")
    foreign_lookup = compute_foreignness_batch(unique_pep_list, proteome_kmers)

    # ── Structural tier 1 ──
    print("  Computing structural tier 1...")
    all_peps_with_wt = []
    all_wts = []
    for windows in mutation_windows.values():
        for pep, wt in windows:
            all_peps_with_wt.append(pep)
            all_wts.append(wt)
    struct_lookup = compute_structural_batch(all_peps_with_wt, all_wts)

    # ── Aggregate: best score per mutation across all windows ──
    print("\n  Aggregating best scores per mutation...")
    best_bigmhc = []
    best_foreign = []
    best_struct = []
    best_pep_per_mut = []

    for idx, row in df.iterrows():
        windows = mutation_windows.get(idx, [])
        if not windows:
            best_bigmhc.append(0.0)
            best_foreign.append(0.5)
            best_struct.append(0.5)
            best_pep_per_mut.append('')
            continue

        # Find the window with the highest BigMHC score
        best_score = -1
        best_pep = ''
        best_f = 0.5
        best_s = 0.5
        for pep, wt in windows:
            score = bigmhc_lookup.get(pep, 0.0)
            if score > best_score:
                best_score = score
                best_pep = pep
                best_f = foreign_lookup.get(pep, 0.5)
                best_s = struct_lookup.get(pep, 0.5)

        best_bigmhc.append(best_score)
        best_foreign.append(best_f)
        best_struct.append(best_s)
        best_pep_per_mut.append(best_pep)

    df['bigmhc_best'] = best_bigmhc
    df['foreign_best'] = best_foreign
    df['struct_best'] = best_struct
    df['best_pep'] = best_pep_per_mut

    # ── Compute other signals ──
    # Agretopicity from pre-computed binding ranks
    agreto = []
    for _, row in df.iterrows():
        mut_rank = row.get('MIN_MUT_RANK_CI_PRIME', row.get('MIN_MUT_RANK_CI_MIXMHC', 50))
        wt_rank = row.get('WT_BEST_RANK_CI_PRIME', row.get('WT_BEST_RANK_CI_MIXMHC', 50))
        try:
            agreto.append(math.log2(max(float(wt_rank), 0.001) / max(float(mut_rank), 0.001)))
        except (ValueError, TypeError):
            agreto.append(0.0)
    df['agretopicity'] = agreto
    df['agreto_norm'] = df['agretopicity'].apply(normalize_agretopicity)
    df['expr_norm'] = df['rnaseq_TPM'].apply(normalize_expression)
    df['ccf_val'] = df['CCF'].fillna(0.5).clip(0, 1)
    df['binding_score'] = 1.0 - df['MIN_MUT_RANK_CI_PRIME'].clip(0, 100) / 100.0

    labels = df['label'].values

    # ══════════════════════════════════════════════
    # BENCHMARK
    # ══════════════════════════════════════════════
    print(f"\n{'='*70}")
    print(f"  MULLER BENCHMARK — {len(df)} mutations ({labels.sum()} immunogenic)")
    print(f"  {df['patient'].nunique()} patients | all signals real | all windows scored")
    print(f"{'='*70}")

    print_metrics("PRIME binding rank (baseline)", labels, df['binding_score'].values)
    print_metrics("BigMHC-IM best window (ours)", labels, df['bigmhc_best'].values)
    print_metrics("Foreignness (k-mer)", labels, df['foreign_best'].values)
    print_metrics("Agretopicity", labels, df['agreto_norm'].values)
    print_metrics("Expression (log TPM)", labels, df['expr_norm'].values)
    print_metrics("Structural tier 1", labels, df['struct_best'].values)
    print_metrics("CCF", labels, df['ccf_val'].values)

    # Composite with all signals
    composite = (
        0.30 * df['bigmhc_best'].values +
        0.15 * df['foreign_best'].values +
        0.10 * df['agreto_norm'].values +
        0.15 * df['binding_score'].values +
        0.10 * df['expr_norm'].values +
        0.10 * df['struct_best'].values +
        0.10 * df['ccf_val'].values
    )
    print_metrics("Our composite (ALL signals, best window)", labels, composite)

    # ── Ablation ──
    print(f"\n{'='*70}")
    print(f"  ABLATION STUDY")
    print(f"{'='*70}")
    ablations = {
        'Full pipeline':         composite,
        '- BigMHC-IM':           composite - 0.30 * df['bigmhc_best'].values,
        '- foreignness':         composite - 0.15 * df['foreign_best'].values,
        '- agretopicity':        composite - 0.10 * df['agreto_norm'].values,
        '- binding':             composite - 0.15 * df['binding_score'].values,
        '- expression':          composite - 0.10 * df['expr_norm'].values,
        '- structural':          composite - 0.10 * df['struct_best'].values,
        '- CCF':                 composite - 0.10 * df['ccf_val'].values,
    }
    print(f"\n  {'Method':<25s} {'AUC':>6s} {'PPV@20':>7s} {'PPV@50':>7s} {'R@20':>6s} {'R@50':>6s}")
    print(f"  {'-'*55}")
    for name, scores in ablations.items():
        a = auc_score(labels, scores)
        print(f"  {name:<25s} {a:>6.3f} {ppv_at_k(labels, scores, 20):>7.3f} "
              f"{ppv_at_k(labels, scores, 50):>7.3f} {recall_at_k(labels, scores, 20):>6.3f} "
              f"{recall_at_k(labels, scores, 50):>6.3f}")

    # ── Per-patient recall ──
    print(f"\n{'='*70}")
    print(f"  PER-PATIENT RECALL")
    print(f"{'='*70}")
    patient_results = []
    for pid in df['patient'].unique():
        pat = df[df['patient'] == pid]
        pat_labels = pat['label'].values
        if pat_labels.sum() == 0:
            continue
        pat_comp = (
            0.30 * pat['bigmhc_best'].values +
            0.15 * pat['foreign_best'].values +
            0.10 * pat['agreto_norm'].values +
            0.15 * pat['binding_score'].values +
            0.10 * pat['expr_norm'].values +
            0.10 * pat['struct_best'].values +
            0.10 * pat['ccf_val'].values
        )
        patient_results.append({
            'patient': pid,
            'n_mut': len(pat),
            'n_pos': int(pat_labels.sum()),
            'R@20': recall_at_k(pat_labels, pat_comp, 20),
            'R@50': recall_at_k(pat_labels, pat_comp, 50),
        })

    pr = pd.DataFrame(patient_results).sort_values('R@20', ascending=False)
    print(f"\n{pr.to_string(index=False)}")
    print(f"\n  Macro Recall@20: {pr['R@20'].mean():.3f}")
    print(f"  Macro Recall@50: {pr['R@50'].mean():.3f}")
    print(f"  Patients R@20 > 0: {(pr['R@20'] > 0).sum()}/{len(pr)}")


if __name__ == '__main__':
    main()
