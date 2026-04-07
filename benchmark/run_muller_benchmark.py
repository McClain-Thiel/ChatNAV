#!/usr/bin/env python3
"""
Benchmark against Muller/Gfeller harmonized dataset (Immunity 2023).

131 patients, 48,306 mutations, 13,483 screened (213 immunogenic).
All signals available: expression (TPM), CCF, clonality, wildtype sequences,
pre-computed binding ranks (NetMHCpan, PRIME, MixMHCpred).

Generates ALL 8-11mer peptide windows per mutation, scores them all with
BigMHC-IM, and takes the best score per mutation — matching how the real
pipeline works.

Outputs the standard artifact set required by AGENT.md §4.
"""

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'bin'))

AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


# ── Data loading ────────────────────────────────────────────────────

def load_muller(split=None, max_patients=None):
    """Load Muller dataset, optionally filtering to a split subset."""
    df = pd.read_csv(PROJECT_ROOT / 'benchmark' / 'muller' / 'Mutation_data_org.txt', sep='\t')
    screened = df[df['response_type'].isin(['CD8', 'negative'])].copy()
    screened['label'] = (screened['response_type'] == 'CD8').astype(int)

    if split:
        splits_path = PROJECT_ROOT / 'benchmark' / 'muller' / 'splits.json'
        if not splits_path.exists():
            raise FileNotFoundError(
                f"splits.json not found at {splits_path}. "
                "Run benchmark/make_splits.py first."
            )
        with open(splits_path) as f:
            splits = json.load(f)
        if split not in ('dev', 'muller_val'):
            raise ValueError(f"Unknown split '{split}'. Use 'dev' or 'muller_val'.")
        patient_ids = splits[split]
        screened = screened[screened['patient'].isin(patient_ids)]

    if max_patients:
        patients = sorted(screened['patient'].unique())[:max_patients]
        screened = screened[screened['patient'].isin(patients)]

    return screened


# ── Peptide window generation ───────────────────────────────────────

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


# ── Scoring helpers ─────────────────────────────────────────────────

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


# ── Metrics ─────────────────────────────────────────────────────────

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


def per_patient_recall_at_k(df, score_col, k):
    """Compute Recall@k per patient, return dict {patient_id: recall}."""
    results = {}
    for pid in df['patient'].unique():
        pat = df[df['patient'] == pid]
        pat_labels = pat['label'].values
        if pat_labels.sum() == 0:
            continue
        pat_scores = pat[score_col].values
        results[pid] = recall_at_k(pat_labels, pat_scores, k)
    return results


def macro_recall_at_k(df, score_col, k):
    """Macro Recall@k: average of per-patient Recall@k."""
    per_patient = per_patient_recall_at_k(df, score_col, k)
    if not per_patient:
        return 0.0
    return np.mean(list(per_patient.values()))


def bootstrap_macro_recall(df, score_col, k, n_bootstrap, seed):
    """Bootstrap CI for macro Recall@k by resampling PATIENTS."""
    rng = np.random.RandomState(seed)

    # Pre-compute per-patient recall
    per_patient = per_patient_recall_at_k(df, score_col, k)
    if not per_patient:
        return 0.0, 0.0, 0.0

    patient_ids = list(per_patient.keys())
    recalls = np.array([per_patient[pid] for pid in patient_ids])
    point_estimate = float(np.mean(recalls))

    boot_means = []
    for _ in range(n_bootstrap):
        sample_idx = rng.choice(len(recalls), size=len(recalls), replace=True)
        boot_means.append(np.mean(recalls[sample_idx]))

    boot_means = np.array(boot_means)
    ci_lower = float(np.percentile(boot_means, 2.5))
    ci_upper = float(np.percentile(boot_means, 97.5))

    return point_estimate, ci_lower, ci_upper


# ── Funnel counting ────────────────────────────────────────────────

def compute_funnel(df, config):
    """Count how many mutations (and how many immunogenic) survive each filter."""
    filters = config.get('hard_filters', {})
    funnel = []

    current = df.copy()
    funnel.append({
        'step': 'all_screened',
        'n_mutations': len(current),
        'n_immunogenic': int(current['label'].sum()),
    })

    # Binding filter
    if 'mhc_binding_rank' in filters:
        thresh = filters['mhc_binding_rank']
        mask = current['MIN_MUT_RANK_CI_PRIME'].fillna(100) / 100.0 <= thresh
        current = current[mask]
        funnel.append({
            'step': f'binding_rank_<={thresh}',
            'n_mutations': len(current),
            'n_immunogenic': int(current['label'].sum()),
        })

    # Expression filter
    if 'min_tpm' in filters and filters['min_tpm'] > 0:
        thresh = filters['min_tpm']
        mask = current['rnaseq_TPM'].fillna(0) >= thresh
        current = current[mask]
        funnel.append({
            'step': f'tpm_>={thresh}',
            'n_mutations': len(current),
            'n_immunogenic': int(current['label'].sum()),
        })

    # CCF filter
    if 'min_ccf' in filters and filters['min_ccf'] > 0:
        thresh = filters['min_ccf']
        mask = current['CCF'].fillna(0) >= thresh
        current = current[mask]
        funnel.append({
            'step': f'ccf_>={thresh}',
            'n_mutations': len(current),
            'n_immunogenic': int(current['label'].sum()),
        })

    return funnel


def apply_hard_filters(df, config):
    """Apply hard filters by setting composite = -inf for filtered mutations.

    Filtered mutations sink below top-K but remain in the denominator for recall,
    so overly aggressive filtering is properly penalized.
    """
    filters = config.get('hard_filters', {})
    mask = np.ones(len(df), dtype=bool)

    if 'mhc_binding_rank' in filters:
        thresh = filters['mhc_binding_rank']
        mask &= (df['MIN_MUT_RANK_CI_PRIME'].fillna(100) / 100.0 <= thresh).values

    if 'min_tpm' in filters and filters['min_tpm'] > 0:
        thresh = filters['min_tpm']
        mask &= (df['rnaseq_TPM'].fillna(0) >= thresh).values

    if 'min_ccf' in filters and filters['min_ccf'] > 0:
        thresh = filters['min_ccf']
        mask &= (df['CCF'].fillna(0) >= thresh).values

    df = df.copy()
    df.loc[~mask, 'composite'] = -np.inf
    return df


# ── Main pipeline ───────────────────────────────────────────────────

def score_mutations(df):
    """Run all scoring on the dataframe. Returns df with score columns added."""
    timings = {}

    # Generate ALL peptide windows per mutation
    t0 = time.time()
    print("\n  Generating all 8-11mer peptide windows per mutation...")
    mutation_windows = {}
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
    timings['window_generation'] = time.time() - t0

    # BigMHC-IM
    t0 = time.time()
    unique_pep_list = list(all_unique_peps)
    hla_list = ['HLA-A*02:01'] * len(unique_pep_list)
    print(f"\n  Running BigMHC-IM on {len(unique_pep_list)} unique peptides...")
    bigmhc_lookup = run_bigmhc_batch(unique_pep_list, hla_list)
    print(f"  BigMHC scored {len(bigmhc_lookup)} peptides")
    timings['bigmhc'] = time.time() - t0

    # Foreignness
    t0 = time.time()
    from score_immunogenicity import build_proteome_kmers
    proteome_path = PROJECT_ROOT / 'reference' / 'proteome' / 'human_proteome.fasta'
    print("  Building proteome k-mer index...")
    proteome_kmers = build_proteome_kmers(str(proteome_path))
    print("  Computing foreignness...")
    foreign_lookup = compute_foreignness_batch(unique_pep_list, proteome_kmers)
    timings['foreignness'] = time.time() - t0

    # Structural tier 1
    t0 = time.time()
    print("  Computing structural tier 1...")
    all_peps_with_wt = []
    all_wts = []
    for windows in mutation_windows.values():
        for pep, wt in windows:
            all_peps_with_wt.append(pep)
            all_wts.append(wt)
    struct_lookup = compute_structural_batch(all_peps_with_wt, all_wts)
    timings['structural'] = time.time() - t0

    # Aggregate: best score per mutation across all windows
    t0 = time.time()
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

    # Other signals
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

    # Differential agretopicity: mutant binds well (< 0.5%) but wildtype doesn't (> 2%)
    mut_rank_col = 'MIN_MUT_RANK_CI_PRIME' if 'MIN_MUT_RANK_CI_PRIME' in df.columns else 'MIN_MUT_RANK_CI_MIXMHC'
    wt_rank_col = 'WT_BEST_RANK_CI_PRIME' if 'WT_BEST_RANK_CI_PRIME' in df.columns else 'WT_BEST_RANK_CI_MIXMHC'
    mut_ranks = pd.to_numeric(df[mut_rank_col], errors='coerce').fillna(50)
    wt_ranks = pd.to_numeric(df[wt_rank_col], errors='coerce').fillna(50)
    df['diff_agreto'] = ((mut_ranks < 0.5) & (wt_ranks > 2.0)).astype(float)

    timings['aggregation'] = time.time() - t0

    return df, timings


DEFAULT_COMPOSITE_WEIGHTS = {
    'bigmhc_best': 0.30,
    'foreign_best': 0.15,
    'agreto_norm': 0.10,
    'binding_score': 0.15,
    'expr_norm': 0.10,
    'struct_best': 0.10,
    'ccf_val': 0.10,
}


def compute_composite(df, config):
    """Compute composite score from config-specified weights."""
    weights = config.get('composite_weights', DEFAULT_COMPOSITE_WEIGHTS)
    composite = np.zeros(len(df))
    for col, w in weights.items():
        if col in df.columns and w != 0:
            composite += w * df[col].values
    df['composite'] = composite
    return df


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


def write_artifacts(df, config, timings, output_dir, n_bootstrap, seed):
    """Write the standard artifact set to output_dir."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    score_col = 'composite'
    labels = df['label'].values

    # recall_at_20.json
    point, ci_lo, ci_hi = bootstrap_macro_recall(df, score_col, 20, n_bootstrap, seed)
    with open(output_dir / 'recall_at_20.json', 'w') as f:
        json.dump({
            'metric': 'macro_recall_at_20',
            'point_estimate': round(point, 4),
            'ci_lower': round(ci_lo, 4),
            'ci_upper': round(ci_hi, 4),
            'n_bootstrap': n_bootstrap,
            'bootstrap_seed': seed,
            'n_patients': len(per_patient_recall_at_k(df, score_col, 20)),
        }, f, indent=2)

    # recall_at_50.json
    point50, ci_lo50, ci_hi50 = bootstrap_macro_recall(df, score_col, 50, n_bootstrap, seed)
    with open(output_dir / 'recall_at_50.json', 'w') as f:
        json.dump({
            'metric': 'macro_recall_at_50',
            'point_estimate': round(point50, 4),
            'ci_lower': round(ci_lo50, 4),
            'ci_upper': round(ci_hi50, 4),
            'n_bootstrap': n_bootstrap,
            'bootstrap_seed': seed,
            'n_patients': len(per_patient_recall_at_k(df, score_col, 50)),
        }, f, indent=2)

    # funnel.json
    funnel = compute_funnel(df, config)
    with open(output_dir / 'funnel.json', 'w') as f:
        json.dump(funnel, f, indent=2)

    # per_patient.tsv
    rows = []
    r20 = per_patient_recall_at_k(df, score_col, 20)
    r50 = per_patient_recall_at_k(df, score_col, 50)
    for pid in sorted(set(list(r20.keys()) + list(r50.keys()))):
        pat = df[df['patient'] == pid]
        rows.append({
            'patient': pid,
            'n_mutations': len(pat),
            'n_immunogenic': int(pat['label'].sum()),
            'recall_at_20': round(r20.get(pid, float('nan')), 4),
            'recall_at_50': round(r50.get(pid, float('nan')), 4),
        })
    pd.DataFrame(rows).to_csv(output_dir / 'per_patient.tsv', sep='\t', index=False)

    # feature_aucs.json
    feature_cols = {
        'bigmhc_best': 'BigMHC-IM',
        'foreign_best': 'Foreignness',
        'agreto_norm': 'Agretopicity',
        'binding_score': 'PRIME binding',
        'expr_norm': 'Expression',
        'struct_best': 'Structural T1',
        'ccf_val': 'CCF',
        'diff_agreto': 'Diff Agretopicity',
        'composite': 'Composite',
    }
    feature_aucs = {}
    for col, name in feature_cols.items():
        if col in df.columns:
            feature_aucs[name] = round(auc_score(labels, df[col].values), 4)
    with open(output_dir / 'feature_aucs.json', 'w') as f:
        json.dump(feature_aucs, f, indent=2)

    # config_snapshot.yaml
    with open(output_dir / 'config_snapshot.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    # git_sha.txt
    try:
        sha = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'],
            cwd=str(PROJECT_ROOT),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        sha = 'unknown'
    with open(output_dir / 'git_sha.txt', 'w') as f:
        f.write(sha + '\n')

    # runtime.json
    with open(output_dir / 'runtime.json', 'w') as f:
        json.dump({k: round(v, 2) for k, v in timings.items()}, f, indent=2)

    return point, ci_lo, ci_hi


def main():
    parser = argparse.ArgumentParser(
        description="Run Muller benchmark with standard artifact output"
    )
    parser.add_argument(
        '--split', choices=['dev', 'muller_val'],
        help='Which patient split to evaluate on (requires splits.json)'
    )
    parser.add_argument(
        '--bootstrap', type=int, default=1000,
        help='Number of bootstrap resamples for CIs (default: 1000)'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for bootstrap (default: 42)'
    )
    parser.add_argument(
        '--max-patients', type=int, default=None,
        help='Limit to first N patients (for quick tests)'
    )
    parser.add_argument(
        '--quick', action='store_true',
        help='Quick mode: 5 patients, 100 bootstrap samples'
    )
    parser.add_argument(
        '--output', type=str, default=None,
        help='Output directory for artifacts (default: print to stdout only)'
    )
    parser.add_argument(
        '--config', type=str,
        default=str(PROJECT_ROOT / 'conf' / 'scoring_weights.yaml'),
        help='Path to scoring weights config'
    )
    parser.add_argument(
        '--profile', type=str, default='high_tmb',
        help='Config profile to use (default: high_tmb)'
    )
    parser.add_argument(
        '--save-scores', type=str, default=None,
        help='Save scored DataFrame to this path (pickle) for reuse'
    )
    parser.add_argument(
        '--load-scores', type=str, default=None,
        help='Load pre-scored DataFrame from this path instead of re-scoring'
    )
    args = parser.parse_args()

    if args.quick:
        args.max_patients = args.max_patients or 5
        args.bootstrap = min(args.bootstrap, 100)

    # Load config
    with open(args.config) as f:
        all_config = yaml.safe_load(f)
    config = all_config.get(args.profile, all_config.get('high_tmb', {}))

    # Load data
    split_label = args.split or 'all'
    print(f"Loading Muller dataset (split={split_label})...")
    df = load_muller(split=args.split, max_patients=args.max_patients)
    print(f"  {len(df)} screened mutations, {df['label'].sum()} immunogenic, "
          f"{df['patient'].nunique()} patients")

    # Score (or load cached scores)
    if args.load_scores:
        print(f"  Loading pre-scored data from {args.load_scores}...")
        df = pd.read_pickle(args.load_scores)
        # Re-apply split/patient filters to cached data
        if args.split:
            splits_path = PROJECT_ROOT / 'benchmark' / 'muller' / 'splits.json'
            with open(splits_path) as f:
                splits = json.load(f)
            df = df[df['patient'].isin(splits[args.split])]
        if args.max_patients:
            patients = sorted(df['patient'].unique())[:args.max_patients]
            df = df[df['patient'].isin(patients)]
        # Recompute derived columns if missing from cache
        if 'diff_agreto' not in df.columns:
            mut_rank_col = 'MIN_MUT_RANK_CI_PRIME' if 'MIN_MUT_RANK_CI_PRIME' in df.columns else 'MIN_MUT_RANK_CI_MIXMHC'
            wt_rank_col = 'WT_BEST_RANK_CI_PRIME' if 'WT_BEST_RANK_CI_PRIME' in df.columns else 'WT_BEST_RANK_CI_MIXMHC'
            mut_ranks = pd.to_numeric(df[mut_rank_col], errors='coerce').fillna(50)
            wt_ranks = pd.to_numeric(df[wt_rank_col], errors='coerce').fillna(50)
            df['diff_agreto'] = ((mut_ranks < 0.5) & (wt_ranks > 2.0)).astype(float)
        timings = {}
        timings['total'] = 0.0
    else:
        total_t0 = time.time()
        df, timings = score_mutations(df)
        timings['total'] = time.time() - total_t0

    if args.save_scores:
        df.to_pickle(args.save_scores)
        print(f"  Saved scored data to {args.save_scores}")

    # Compute composite from config weights
    df = compute_composite(df, config)

    # Apply hard filters: filtered mutations get composite=-inf so they
    # sink below top-K, but remain in the denominator for recall
    df = apply_hard_filters(df, config)
    n_surviving = (df['composite'] > -np.inf).sum()
    n_pos_surviving = df.loc[df['composite'] > -np.inf, 'label'].sum()
    print(f"\n  After hard filters: {n_surviving}/{len(df)} mutations survive "
          f"({n_pos_surviving}/{df['label'].sum()} immunogenic)")

    labels = df['label'].values

    # Print report
    print(f"\n{'='*70}")
    print(f"  MULLER BENCHMARK — {len(df)} mutations ({labels.sum()} immunogenic)")
    print(f"  {df['patient'].nunique()} patients | split={split_label}")
    print(f"{'='*70}")

    print_metrics("PRIME binding rank (baseline)", labels, df['binding_score'].values)
    print_metrics("BigMHC-IM best window (ours)", labels, df['bigmhc_best'].values)
    print_metrics("Foreignness (k-mer)", labels, df['foreign_best'].values)
    print_metrics("Agretopicity", labels, df['agreto_norm'].values)
    print_metrics("Expression (log TPM)", labels, df['expr_norm'].values)
    print_metrics("Structural tier 1", labels, df['struct_best'].values)
    print_metrics("CCF", labels, df['ccf_val'].values)
    if 'diff_agreto' in df.columns:
        print_metrics("Diff agretopicity", labels, df['diff_agreto'].values)
    print_metrics("Our composite (ALL signals, best window)", labels, df['composite'].values)

    # Bootstrap macro Recall@20
    print(f"\n{'='*70}")
    print(f"  MACRO RECALL (bootstrap={args.bootstrap}, seed={args.seed})")
    print(f"{'='*70}")

    point20, lo20, hi20 = bootstrap_macro_recall(
        df, 'composite', 20, args.bootstrap, args.seed
    )
    point50, lo50, hi50 = bootstrap_macro_recall(
        df, 'composite', 50, args.bootstrap, args.seed
    )
    print(f"\n  Macro Recall@20: {point20:.4f} [{lo20:.4f}, {hi20:.4f}]")
    print(f"  Macro Recall@50: {point50:.4f} [{lo50:.4f}, {hi50:.4f}]")

    # Per-patient summary
    print(f"\n{'='*70}")
    print(f"  PER-PATIENT RECALL")
    print(f"{'='*70}")
    r20 = per_patient_recall_at_k(df, 'composite', 20)
    r50 = per_patient_recall_at_k(df, 'composite', 50)
    patient_rows = []
    for pid in sorted(r20.keys()):
        pat = df[df['patient'] == pid]
        patient_rows.append({
            'patient': pid,
            'n_mut': len(pat),
            'n_pos': int(pat['label'].sum()),
            'R@20': r20.get(pid, 0.0),
            'R@50': r50.get(pid, 0.0),
        })
    pr = pd.DataFrame(patient_rows).sort_values('R@20', ascending=False)
    print(f"\n{pr.to_string(index=False)}")
    print(f"\n  Patients with R@20 > 0: {(pr['R@20'] > 0).sum()}/{len(pr)}")

    # Write artifacts
    if args.output:
        print(f"\n  Writing artifacts to {args.output}...")
        write_artifacts(df, config, timings, args.output, args.bootstrap, args.seed)
        print(f"  Done. Artifacts written to {args.output}")

    return point20


if __name__ == '__main__':
    main()
