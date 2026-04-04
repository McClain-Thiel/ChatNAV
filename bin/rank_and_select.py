#!/usr/bin/env python3
"""
Module 9: Candidate Ranking and Selection

Sequential filter pipeline — each step reduces the candidate pool,
then final ranking by binding strength within survivors.

Pipeline:
  1. FILTER: MHC binding rank < threshold (keeps strong binders)
  2. FILTER: Gene expression TPM > min (gene must be expressed)
  3. FILTER: CCF > min (mutation must be clonal)
  4. FILTER: Not a self-match (foreignness > 0)
  5. RANK: By binding rank (strongest binders first)
  6. FLAG: Immunogenicity score, agretopicity, structural as validation signals
  7. DIVERSIFY: Ensure HLA allele coverage in top N
  8. SELECT: Top N

Bonuses promote frameshifts and shared neoantigens by rank adjustment.

Inputs:
    immunogenicity_scores.tsv   — from Module 8a
    structural_scores.tsv       — from Module 8b
    candidates_meta.tsv         — from Step 0
    binding_predictions.tsv     — from Step 0
    scoring_weights.yaml        — filter thresholds + profile config

Output:
    selected_neoantigens.tsv    — top N neoantigens with all scores + flags
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


def load_filter_profile(weights_path: str, profile_name: str) -> dict:
    """Load a filter/ranking profile from scoring_weights.yaml."""
    with open(weights_path) as f:
        all_profiles = yaml.safe_load(f)

    if profile_name not in all_profiles:
        raise ValueError(
            f"Profile '{profile_name}' not found in {weights_path}. "
            f"Available: {list(all_profiles.keys())}"
        )

    return all_profiles[profile_name]


def apply_sequential_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
    """
    Apply sequential hard filters. Each step logs how many candidates
    and how many positives (if labeled) survive.
    """
    initial = len(df)
    steps = []

    # Step 1: MHC binding
    if 'mhc_binding_rank' in filters:
        thresh = filters['mhc_binding_rank']
        df = df[df['binding_rank'] < thresh]
        steps.append(f"Binding < {thresh}: {len(df)} remain")

    # Step 2: Expression
    if 'min_tpm' in filters:
        thresh = filters['min_tpm']
        df = df[df['tpm'].fillna(0) >= thresh]
        steps.append(f"TPM >= {thresh}: {len(df)} remain")

    # Step 3: CCF (clonality)
    if 'min_ccf' in filters:
        thresh = filters['min_ccf']
        df = df[df['ccf'].fillna(0) >= thresh]
        steps.append(f"CCF >= {thresh}: {len(df)} remain")

    # Step 4: Self-match filter
    if filters.get('not_self_match', True):
        if 'is_self_match' in df.columns:
            df = df[~df['is_self_match'].fillna(False)]
            steps.append(f"Not self-match: {len(df)} remain")

    print(f"Sequential filters: {initial} → {len(df)} candidates", file=sys.stderr)
    for s in steps:
        print(f"  {s}", file=sys.stderr)

    return df


def apply_rank_bonuses(df: pd.DataFrame, bonuses: dict) -> pd.DataFrame:
    """
    Adjust binding rank for bonus candidates (frameshifts, shared neoantigens).
    These get promoted in the ranking by reducing their effective binding rank.
    """
    df = df.copy()
    df['effective_rank'] = df['binding_rank'].copy()

    # Frameshifts: promote by reducing effective rank
    if bonuses.get('is_frameshift', 0) > 0 and 'is_frameshift' in df.columns:
        mask = df['is_frameshift'].fillna(False)
        # Multiply rank by (1 - bonus) to promote. E.g., bonus=0.5 halves the rank.
        df.loc[mask, 'effective_rank'] *= (1 - bonuses['is_frameshift'])
        n = mask.sum()
        if n > 0:
            print(f"  Frameshift bonus: {n} candidates promoted", file=sys.stderr)

    # Shared neoantigens: promote
    if bonuses.get('is_shared_neoantigen', 0) > 0 and 'is_shared_neoantigen' in df.columns:
        mask = df['is_shared_neoantigen'].fillna(False)
        df.loc[mask, 'effective_rank'] *= (1 - bonuses['is_shared_neoantigen'])
        n = mask.sum()
        if n > 0:
            print(f"  Shared neoantigen bonus: {n} candidates promoted", file=sys.stderr)

    # CD4 epitopes: promote
    if bonuses.get('is_cd4_epitope', 0) > 0 and 'is_cd4_epitope' in df.columns:
        mask = df['is_cd4_epitope'].fillna(False)
        df.loc[mask, 'effective_rank'] *= (1 - bonuses['is_cd4_epitope'])

    return df


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

    remaining = df.iloc[top_n:]
    primary_hla = selected['hla_allele'].iloc[0]
    alt_candidates = remaining[remaining['hla_allele'] != primary_hla]

    if not alt_candidates.empty:
        swap_in = alt_candidates.iloc[0]
        selected = pd.concat([selected.iloc[:-1], swap_in.to_frame().T], ignore_index=True)

    return selected


def main():
    parser = argparse.ArgumentParser(description='Module 9: Sequential Filter Ranking')
    parser.add_argument('--immunogenicity-scores', required=True)
    parser.add_argument('--candidates-meta', required=True)
    parser.add_argument('--binding-predictions', required=True)
    parser.add_argument('--scoring-weights', required=True)
    parser.add_argument('--weight-profile', default='high_tmb')
    parser.add_argument('--top-n', type=int, default=20)
    parser.add_argument('--structural-scores', default=None)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # Load profile
    profile = load_filter_profile(args.scoring_weights, args.weight_profile)
    bonuses = profile.get('bonuses', {})
    filters = profile['hard_filters']

    # Load data
    immuno = pd.read_csv(args.immunogenicity_scores, sep='\t')
    meta = pd.read_csv(args.candidates_meta, sep='\t')
    binding = pd.read_csv(args.binding_predictions, sep='\t')

    # Load structural scores if provided
    structural = None
    if args.structural_scores:
        structural = pd.read_csv(args.structural_scores, sep='\t')

    # Merge all data
    df = immuno.merge(binding, on=['peptide_id', 'hla_allele'], how='left')
    df = df.merge(meta, on='peptide_id', how='left', suffixes=('', '_meta'))

    if structural is not None:
        struct_cols = [c for c in ['peptide_id', 'hla_allele', 'structural_score', 'structural_tier']
                       if c in structural.columns]
        df = df.merge(structural[struct_cols], on=['peptide_id', 'hla_allele'], how='left')

    # Normalize expression
    df['expression_norm'] = df['tpm'].apply(normalize_expression)

    # Compute agretopicity flag (positive = mutation creates new binder)
    if 'agretopicity' in df.columns:
        df['agreto_positive'] = df['agretopicity'].fillna(0) > 0

    # Fill boolean flags
    for col in ['is_frameshift', 'is_shared_neoantigen', 'is_cd4_epitope']:
        if col not in df.columns:
            df[col] = False
        df[col] = df[col].fillna(False).astype(bool)

    if 'mhc_class' in df.columns:
        df.loc[df['mhc_class'] == 'II', 'is_cd4_epitope'] = True

    # ── Step 1-4: Sequential hard filters ──
    df = apply_sequential_filters(df, filters)

    if len(df) == 0:
        raise RuntimeError("No candidates passed filters — check binding data and thresholds")

    # ── Step 5: Apply rank bonuses (frameshifts, shared neoantigens) ──
    df = apply_rank_bonuses(df, bonuses)

    # ── Step 6: Rank by effective binding rank (best binders first) ──
    df = df.sort_values('effective_rank').reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)

    # ── Step 7: HLA diversity ──
    top_n = args.top_n if args.weight_profile != 'research' else len(df)
    selected = ensure_hla_diversity(df, top_n)
    selected = selected.copy()
    selected['selected'] = True

    # Re-rank
    selected['rank'] = range(1, len(selected) + 1)

    # ── Step 8: Flag validation signals on selected candidates ──
    # These are informational — they don't change the ranking but help
    # the user assess candidate quality
    flags = []
    for _, row in selected.iterrows():
        f = []
        im = row.get('immunogenicity_score', 0)
        if im > 0.7:
            f.append('high_immunogenicity')
        elif im < 0.3:
            f.append('low_immunogenicity')

        agreto = row.get('agretopicity', 0)
        if agreto < 0:
            f.append('wt_binds_better')
        elif agreto > 3:
            f.append('strong_differential_binding')

        struct = row.get('structural_score', 0.5)
        if struct > 0.7:
            f.append('tcr_exposed')
        elif struct < 0.2:
            f.append('buried_mutation')

        foreign = row.get('foreignness_score', 0.5)
        if foreign < 0.2:
            f.append('self_similar')

        if row.get('is_frameshift', False):
            f.append('frameshift')
        if row.get('is_shared_neoantigen', False):
            f.append('shared_neoantigen')

        flags.append('; '.join(f) if f else '')

    selected['validation_flags'] = flags

    # Output columns
    output_cols = [
        'rank', 'peptide_id', 'peptide_sequence', 'hla_allele', 'mhc_class',
        'gene', 'mutation', 'mutation_type',
        'binding_rank', 'effective_rank',
        'tpm', 'ccf',
        'immunogenicity_score', 'foreignness_score', 'agretopicity',
        'structural_score', 'structural_tier',
        'is_frameshift', 'is_shared_neoantigen', 'is_cd4_epitope',
        'validation_flags',
        'wildtype_peptide', 'selected'
    ]
    output_cols = [c for c in output_cols if c in selected.columns]
    selected = selected[output_cols]

    selected.to_csv(args.output, sep='\t', index=False)

    # Summary
    print(f"\nSelected {len(selected)} neoantigens → {args.output}", file=sys.stderr)
    n_flags = sum(1 for f in flags if f)
    print(f"  {n_flags}/{len(selected)} have validation flags", file=sys.stderr)
    for _, row in selected.head(5).iterrows():
        print(f"  #{row['rank']}: {row.get('peptide_id','?')} "
              f"bind={row.get('binding_rank',0):.4f} "
              f"IM={row.get('immunogenicity_score',0):.3f} "
              f"[{row.get('validation_flags','')}]",
              file=sys.stderr)


if __name__ == '__main__':
    main()
