#!/usr/bin/env python3
"""
Module 9: Candidate Ranking and Selection

Integrates all signals into a composite score and selects the top N
neoantigens for the vaccine construct.

Composite score formula (configurable weights):
    composite = (w1 × immunogenicity + w2 × foreignness + w3 × agretopicity
                 + w4 × (1 - binding_rank) + w5 × stability + w6 × norm_expression
                 + w7 × ccf + w8 × structural)
    + bonuses (frameshift, shared neoantigen, CD4+ epitope)

Hard filters applied before scoring:
    binding_rank < threshold, TPM > min, CCF > min, agretopicity > 0

Inputs:
    immunogenicity_scores.tsv   — from Module 8
    candidates_meta.tsv         — from Module 6
    binding_predictions.tsv     — from Module 7
    scoring_weights.yaml        — weight profile config

Output:
    selected_neoantigens.tsv    — top N neoantigens with all scores
"""

import argparse
import math
import sys

import pandas as pd
import yaml


def normalize_expression(tpm: float) -> float:
    """Normalize TPM to [0, 1] range: log10(tpm + 1) / 5"""
    if tpm is None or pd.isna(tpm) or tpm < 0:
        return 0.0
    return min(math.log10(tpm + 1) / 5.0, 1.0)


def normalize_agretopicity(agreto: float) -> float:
    """Normalize agretopicity to [0, 1] range using sigmoid-like mapping."""
    if agreto is None or pd.isna(agreto):
        return 0.0
    # Agretopicity of 5.0 maps to ~0.95; 0.0 maps to 0.5; negative maps lower
    return 1.0 / (1.0 + math.exp(-0.5 * agreto))


def load_weight_profile(weights_path: str, profile_name: str) -> dict:
    """Load a weight profile from the scoring_weights.yaml config."""
    with open(weights_path) as f:
        all_profiles = yaml.safe_load(f)

    if profile_name not in all_profiles:
        raise ValueError(
            f"Weight profile '{profile_name}' not found in {weights_path}. "
            f"Available profiles: {list(all_profiles.keys())}"
        )

    return all_profiles[profile_name]


def apply_hard_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
    """Apply hard filters; return only rows that pass all filters."""
    initial_count = len(df)

    if 'mhc_binding_rank' in filters:
        df = df[df['binding_rank'] < filters['mhc_binding_rank']]

    if 'min_tpm' in filters:
        df = df[df['tpm'].fillna(0) >= filters['min_tpm']]

    if 'min_ccf' in filters:
        df = df[df['ccf'].fillna(0) >= filters['min_ccf']]

    if 'min_agretopicity' in filters:
        df = df[df['agretopicity'].fillna(0) > filters['min_agretopicity']]

    if filters.get('not_self_match', True):
        if 'is_self_match' in df.columns:
            df = df[~df['is_self_match'].fillna(False)]

    print(f"Hard filters: {initial_count} → {len(df)} candidates", file=sys.stderr)
    return df


def compute_composite_scores(df: pd.DataFrame, weights: dict, bonuses: dict,
                              structural_enabled: bool = False) -> pd.DataFrame:
    """Compute composite score for each candidate."""
    # Validate required columns exist and have no missing values
    required_cols = ['immunogenicity_score', 'foreignness_score', 'agretopicity_norm',
                     'binding_rank', 'expression_norm']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' missing from merged data")
        n_missing = df[col].isna().sum()
        if n_missing > 0:
            raise ValueError(f"Column '{col}' has {n_missing} missing values — fix upstream data")

    w = weights.copy()
    if structural_enabled and w.get('structural', 0) == 0:
        w['structural'] = 0.10
        w['immunogenicity_score'] = max(0, w.get('immunogenicity_score', 0.35) - 0.10)

    # structural_score and stability_rank are optional (default 0 / 0.5 if not present)
    stability = (1 - df['stability_rank']) if 'stability_rank' in df.columns and not df['stability_rank'].isna().all() else 0.5
    structural = df['structural_score'].fillna(0) if 'structural_score' in df.columns else 0

    scores = (
        w.get('immunogenicity_score', 0.35) * df['immunogenicity_score']
        + w.get('foreignness_score', 0.15) * df['foreignness_score']
        + w.get('agretopicity', 0.10) * df['agretopicity_norm']
        + w.get('mhc_binding', 0.15) * (1 - df['binding_rank'])
        + w.get('stability', 0.10) * stability
        + w.get('expression', 0.10) * df['expression_norm']
        + w.get('ccf', 0.05) * df['ccf'].fillna(1.0)
        + w.get('structural', 0.00) * structural
    )

    # Add bonuses
    if bonuses.get('is_frameshift', 0) > 0:
        scores += bonuses['is_frameshift'] * df['is_frameshift'].fillna(False).astype(float)
    if bonuses.get('is_shared_neoantigen', 0) > 0:
        scores += bonuses['is_shared_neoantigen'] * df['is_shared_neoantigen'].fillna(False).astype(float)
    if bonuses.get('is_cd4_epitope', 0) > 0:
        scores += bonuses['is_cd4_epitope'] * df['is_cd4_epitope'].fillna(False).astype(float)

    return scores


def ensure_hla_diversity(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """
    Select top N with HLA allele diversity: ensure at least 2 different
    HLA alleles are represented if possible.
    """
    if len(df) <= top_n:
        return df

    selected = df.head(top_n).copy()
    unique_hla = selected['hla_allele'].nunique()

    if unique_hla >= 2:
        return selected

    # Try to swap in a candidate with a different HLA allele
    remaining = df.iloc[top_n:]
    primary_hla = selected['hla_allele'].iloc[0]
    alt_candidates = remaining[remaining['hla_allele'] != primary_hla]

    if not alt_candidates.empty:
        # Swap the lowest-scoring selected candidate with the best alternative
        swap_in = alt_candidates.iloc[0]
        selected = pd.concat([selected.iloc[:-1], swap_in.to_frame().T], ignore_index=True)

    return selected


def main():
    parser = argparse.ArgumentParser(description='Module 9: Candidate Ranking and Selection')
    parser.add_argument('--immunogenicity-scores', required=True)
    parser.add_argument('--candidates-meta', required=True)
    parser.add_argument('--binding-predictions', required=True)
    parser.add_argument('--scoring-weights', required=True)
    parser.add_argument('--weight-profile', default='high_tmb')
    parser.add_argument('--top-n', type=int, default=20)
    parser.add_argument('--structural-scores', default=None,
                        help='structural_scores.tsv from Module 8b')
    parser.add_argument('--output', required=True)
    parser.add_argument('--output-stats', default=None,
                        help='Optional: write pipeline_stats.json with funnel counts')
    args = parser.parse_args()

    # Load weight profile
    profile = load_weight_profile(args.scoring_weights, args.weight_profile)
    weights = profile['weights']
    bonuses = profile['bonuses']
    filters = profile['hard_filters']

    # Load data
    immuno = pd.read_csv(args.immunogenicity_scores, sep='\t')
    meta = pd.read_csv(args.candidates_meta, sep='\t')
    binding = pd.read_csv(args.binding_predictions, sep='\t')

    # Load structural scores if provided
    structural = None
    if args.structural_scores:
        structural = pd.read_csv(args.structural_scores, sep='\t')
        print(f"Loaded {len(structural)} structural scores", file=sys.stderr)

    # Merge all data on peptide_id + hla_allele
    df = immuno.merge(binding, on=['peptide_id', 'hla_allele'], how='left')
    df = df.merge(meta, on='peptide_id', how='left', suffixes=('', '_meta'))

    # Merge structural scores
    if structural is not None:
        struct_cols = ['peptide_id', 'hla_allele', 'structural_score', 'structural_tier']
        struct_cols = [c for c in struct_cols if c in structural.columns]
        df = df.merge(structural[struct_cols], on=['peptide_id', 'hla_allele'], how='left')

    # Track funnel counts
    n_total_scored = len(df)
    n_unique_mutations = df['peptide_id'].apply(lambda x: '_'.join(x.rsplit('_', 1)[:-1]) if '_' in x else x).nunique()

    # Normalize expression and agretopicity
    df['expression_norm'] = df['tpm'].apply(normalize_expression)
    df['agretopicity_norm'] = df['agretopicity'].apply(normalize_agretopicity)

    # Fill missing boolean flags
    for col in ['is_frameshift', 'is_shared_neoantigen', 'is_cd4_epitope']:
        if col not in df.columns:
            df[col] = False
        df[col] = df[col].fillna(False).astype(bool)

    # Detect MHC-II binders as CD4+ epitopes
    if 'mhc_class' in df.columns:
        df.loc[df['mhc_class'] == 'II', 'is_cd4_epitope'] = True

    # Save pre-filter state for stats, count binders
    df_pre_filter = df.copy()

    # Apply hard filters
    df = apply_hard_filters(df, filters)
    n_after_filter = len(df)

    if len(df) == 0:
        raise RuntimeError("No candidates passed hard filters — check binding data and filter thresholds")

    # Compute composite scores
    has_structural = 'structural_score' in df.columns and not df['structural_score'].isna().all()
    df['composite_score'] = compute_composite_scores(
        df, weights, bonuses, structural_enabled=has_structural
    )

    # Rank by composite score
    df = df.sort_values('composite_score', ascending=False).reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)

    # Select top N with HLA diversity
    top_n = args.top_n if args.weight_profile != 'research' else len(df)
    selected = ensure_hla_diversity(df, top_n)
    selected = selected.copy()
    selected['selected'] = True

    # Output columns
    output_cols = [
        'rank', 'peptide_id', 'peptide_sequence', 'hla_allele', 'mhc_class',
        'gene', 'mutation', 'mutation_type',
        'binding_rank', 'stability_rank',
        'tpm', 'expression_imputed', 'ccf',
        'immunogenicity_score', 'foreignness_score', 'agretopicity',
        'structural_score', 'structural_tier', 'composite_score',
        'is_frameshift', 'is_shared_neoantigen', 'is_cd4_epitope',
        'wildtype_peptide', 'selected'
    ]
    # Keep only columns that exist
    output_cols = [c for c in output_cols if c in selected.columns]
    selected = selected[output_cols]

    # Re-rank after selection
    selected['rank'] = range(1, len(selected) + 1)

    selected.to_csv(args.output, sep='\t', index=False)
    print(f"Selected {len(selected)} neoantigens → {args.output}")

    # Write pipeline stats if requested
    if args.output_stats:
        import json as _json
        n_binders = (df_pre_filter['binding_rank'] < 0.02).sum() if 'df_pre_filter' in dir() else len(df)
        stats = {
            'n_windows': n_total_scored,
            'n_unique_mutations': n_unique_mutations,
            'n_binders': int(n_binders),
            'n_filtered': int(n_after_filter),
            'n_selected': len(selected),
        }
        with open(args.output_stats, 'w') as f:
            _json.dump(stats, f, indent=2)
        print(f"Pipeline stats → {args.output_stats}", file=sys.stderr)


if __name__ == '__main__':
    main()
