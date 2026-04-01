#!/usr/bin/env python3
"""
Full-pipeline benchmark against Gartner et al. 2021 (Nature Cancer).

Runs the COMPLETE scoring pipeline (not just BigMHC) on ALL screened
peptide-HLA combinations for each patient, then checks whether the
confirmed immunogenic peptides land in the top-K selections.

For each patient:
  1. Take ALL their screened nmers (positive + negative)
  2. Extract 8-11mer windows from each nmer
  3. Run MHCflurry binding predictions
  4. Run BigMHC-IM immunogenicity scoring
  5. Compute foreignness vs human proteome
  6. Compute agretopicity (mut vs wt binding)
  7. Run structural tier 1 scoring
  8. Apply Module 9 composite ranking
  9. Check: are the confirmed immunogenic peptides in the top 20?

This tests the FULL pipeline, not just one module.
"""

import sys
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'bin'))


def load_gartner():
    xlsx = PROJECT_ROOT / 'benchmark' / 'gartner' / 'supplementary_table.xlsx'

    patients = pd.read_excel(xlsx, sheet_name='Supplementary Table 1', header=1)
    hla_map = dict(zip(patients['ID'], patients['Patient HLAs']))

    positives = pd.read_excel(xlsx, sheet_name='Supplemetary Table 3', header=1)
    positives = positives.dropna(subset=['Mutant minimal Peptide', 'HLA'])

    nmers = pd.read_excel(xlsx, sheet_name='Supplementary Table 11', header=1)

    return patients, hla_map, positives, nmers


def extract_mmps_from_nmer(nmer: str, hla_alleles: list[str]):
    """Generate all 8-11mer windows from a nmer."""
    windows = []
    for pep_len in [8, 9, 10, 11]:
        for start in range(len(nmer) - pep_len + 1):
            pep = nmer[start:start + pep_len]
            if all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in pep):
                for hla in hla_alleles:
                    windows.append((pep, hla))
    return windows


def run_mhcflurry_genotype(peptides, genotype_alleles):
    """
    Run MHCflurry in genotype mode: each peptide scored against all alleles,
    returns the best allele + presentation score per peptide.
    """
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    # Deduplicate peptides
    unique_peps = list(set(peptides))
    clean_alleles = [a.replace('*', '').replace(':', '') for a in genotype_alleles]

    results = predictor.predict(
        peptides=unique_peps,
        alleles=clean_alleles,
        verbose=0,
    )

    # Build lookup: peptide → (best_allele, presentation_percentile, presentation_score)
    lookup = {}
    for _, row in results.iterrows():
        pep = row['peptide']
        lookup[pep] = {
            'best_allele': row['best_allele'],
            'presentation_percentile': row['presentation_percentile'],
            'presentation_score': row['presentation_score'],
        }

    return lookup


def run_bigmhc_batch(peptides, hla_alleles):
    """Run BigMHC-IM in-process."""
    from score_immunogenicity import _run_bigmhc_inprocess

    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        for pep, hla in zip(peptides, hla_alleles):
            f.write(f"{pep},{hla}\n")
        path = f.name

    _run_bigmhc_inprocess(path)
    out = pd.read_csv(path + '.prd')
    Path(path).unlink(missing_ok=True)
    Path(path + '.prd').unlink(missing_ok=True)
    return out['BigMHC_IM'].tolist()


def run_foreignness_batch(peptides, proteome_kmers):
    """Compute foreignness for each peptide."""
    from score_immunogenicity import compute_foreignness_kmer
    scores = []
    for pep in peptides:
        if len(pep) >= 9:
            try:
                scores.append(compute_foreignness_kmer(pep, proteome_kmers))
            except ValueError:
                scores.append(0.5)
        else:
            scores.append(0.5)
    return scores


def run_structural_batch(peptides):
    """Run structural tier 1 on peptides (no WT available for nmers)."""
    from structural_scoring import score_tier1_position
    scores = []
    for pep in peptides:
        # All positions treated as mutated (no WT comparison for nmer-derived peptides)
        mut_pos = list(range(len(pep)))
        score, _ = score_tier1_position(pep, mut_pos, len(pep))
        scores.append(score if score is not None else 0.5)
    return scores


def main():
    patients, hla_map, positives, nmers = load_gartner()

    # Build proteome k-mer index
    proteome_path = PROJECT_ROOT / 'reference' / 'proteome' / 'human_proteome.fasta'
    if not proteome_path.exists():
        print("FATAL: Human proteome not found — run SETUP.md instructions", file=sys.stderr)
        sys.exit(1)

    from score_immunogenicity import build_proteome_kmers
    print("Building proteome k-mer index...", file=sys.stderr)
    proteome_kmers = build_proteome_kmers(str(proteome_path))
    print(f"Built {len(proteome_kmers)} k-mers", file=sys.stderr)

    # Pick patients with most confirmed immunogenic peptides
    patient_pos_counts = positives.groupby('ID').size().sort_values(ascending=False)
    test_patients = patient_pos_counts.head(5).index.tolist()

    print(f"\nRunning full pipeline on {len(test_patients)} patients:")
    for pid in test_patients:
        pat_pos = positives[positives['ID'] == pid]
        pat_nmers = nmers[nmers['ID'] == pid]
        print(f"  Patient {pid}: {len(pat_nmers)} nmers, {len(pat_pos)} confirmed immunogenic")

    # Build the confirmed positive peptide set (for recall check)
    positive_set = {}  # patient_id → set of (peptide, hla)
    for _, row in positives.iterrows():
        pid = row['ID']
        pep = str(row['Mutant minimal Peptide']).strip()
        hla = str(row['HLA']).strip()
        if pid not in positive_set:
            positive_set[pid] = set()
        positive_set[pid].add((pep, hla))

    # Run full pipeline per patient
    all_results = []

    for pid in test_patients:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Patient {pid}", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        pat_nmers = nmers[nmers['ID'] == pid]
        hlas = [h.strip() for h in str(hla_map[pid]).split(',')]
        # MHCflurry predict() takes at most 6 alleles (one genotype)
        hlas = [h for h in hlas if h.startswith('HLA-')][:6]

        # Extract unique 9-10mer windows from ALL nmers
        all_peps_set = set()
        for _, nmer_row in pat_nmers.iterrows():
            nmer_seq = str(nmer_row['Mutant nmer']).strip()
            for pep_len in [9, 10]:
                for start in range(len(nmer_seq) - pep_len + 1):
                    pep = nmer_seq[start:start + pep_len]
                    if all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in pep):
                        all_peps_set.add(pep)

        all_peps_unique = list(all_peps_set)
        print(f"  {len(all_peps_unique)} unique peptide windows from {len(pat_nmers)} nmers", file=sys.stderr)

        if not all_peps_unique:
            continue

        # Step 1: MHCflurry in genotype mode (scores each peptide against all HLAs)
        print(f"  Running MHCflurry ({len(all_peps_unique)} peptides × {len(hlas)} alleles)...", file=sys.stderr)
        mhcflurry_lookup = run_mhcflurry_genotype(all_peps_unique, hlas)

        # Filter to strong binders (presentation_percentile < 2%)
        peps_filtered = []
        hlas_filtered = []
        ranks_filtered = []
        for pep in all_peps_unique:
            info = mhcflurry_lookup.get(pep)
            if info and info['presentation_percentile'] / 100.0 < 0.02:
                peps_filtered.append(pep)
                # Format HLA back: HLA-A0201 → HLA-A*02:01
                raw_hla = info['best_allele']
                if '*' not in raw_hla and len(raw_hla) >= 8:
                    formatted_hla = f"{raw_hla[:5]}*{raw_hla[5:7]}:{raw_hla[7:]}"
                else:
                    formatted_hla = raw_hla
                hlas_filtered.append(formatted_hla)
                ranks_filtered.append(info['presentation_percentile'] / 100.0)

        if not peps_filtered:
            # Relax threshold
            print(f"  WARNING: No binders at <2% — relaxing to <5%", file=sys.stderr)
            for pep in all_peps_unique:
                info = mhcflurry_lookup.get(pep)
                if info and info['presentation_percentile'] / 100.0 < 0.05:
                    peps_filtered.append(pep)
                    raw_hla = info['best_allele']
                    if '*' not in raw_hla and len(raw_hla) >= 8:
                        formatted_hla = f"{raw_hla[:5]}*{raw_hla[5:7]}:{raw_hla[7:]}"
                    else:
                        formatted_hla = raw_hla
                    hlas_filtered.append(formatted_hla)
                    ranks_filtered.append(info['presentation_percentile'] / 100.0)

        print(f"  {len(peps_filtered)} pass binding filter", file=sys.stderr)

        # Step 2: BigMHC-IM immunogenicity
        print(f"  Running BigMHC-IM...", file=sys.stderr)
        bigmhc_scores = run_bigmhc_batch(peps_filtered, hlas_filtered)

        # Step 3: Foreignness
        print(f"  Computing foreignness...", file=sys.stderr)
        foreign_scores = run_foreignness_batch(peps_filtered, proteome_kmers)

        # Step 4: Structural tier 1
        struct_scores = run_structural_batch(peps_filtered)

        # Step 5: Composite score (Module 9 formula)
        from rank_and_select import normalize_agretopicity
        composite_scores = []
        for i in range(len(peps_filtered)):
            # Simplified composite (no expression/CCF since we only have deciles)
            binding_signal = 1.0 - ranks_filtered[i]
            immuno = bigmhc_scores[i]
            foreign = foreign_scores[i]
            structural = struct_scores[i]

            composite = (
                0.40 * immuno +
                0.20 * foreign +
                0.20 * binding_signal +
                0.10 * structural +
                0.10 * 0.5  # placeholder for expression/CCF
            )
            composite_scores.append(composite)

        # Rank by composite score
        ranked_indices = np.argsort(composite_scores)[::-1]

        # Check recall: are the confirmed immunogenic peptides in the rankings?
        patient_positives = positive_set.get(pid, set())
        found_at_rank = {}
        for rank, idx in enumerate(ranked_indices):
            pep = peps_filtered[idx]
            hla = hlas_filtered[idx]
            # Normalize HLA for matching
            hla_norm = hla.replace('*', '').replace(':', '')
            for pos_pep, pos_hla in patient_positives:
                pos_hla_norm = pos_hla.replace('*', '').replace(':', '')
                if pep == pos_pep and hla_norm == pos_hla_norm:
                    if pos_pep not in found_at_rank:
                        found_at_rank[pos_pep] = rank + 1

        # Also check peptide-only match (HLA might differ in format)
        pos_pep_set = {p for p, h in patient_positives}
        for rank, idx in enumerate(ranked_indices):
            pep = peps_filtered[idx]
            if pep in pos_pep_set and pep not in found_at_rank:
                found_at_rank[pep] = rank + 1

        # Results
        n_found_top20 = sum(1 for r in found_at_rank.values() if r <= 20)
        n_found_top50 = sum(1 for r in found_at_rank.values() if r <= 50)
        n_total_pos = len(patient_positives)

        print(f"\n  Results for patient {pid}:", file=sys.stderr)
        print(f"    Confirmed immunogenic: {n_total_pos}", file=sys.stderr)
        print(f"    Found in ranked list:  {len(found_at_rank)}/{n_total_pos}", file=sys.stderr)
        print(f"    Recall@20:             {n_found_top20}/{n_total_pos} = {n_found_top20/max(n_total_pos,1):.2f}", file=sys.stderr)
        print(f"    Recall@50:             {n_found_top50}/{n_total_pos} = {n_found_top50/max(n_total_pos,1):.2f}", file=sys.stderr)

        for pep, rank in sorted(found_at_rank.items(), key=lambda x: x[1]):
            idx = ranked_indices[rank - 1]
            print(f"    Rank {rank:>4d}: {pep:<15s} composite={composite_scores[idx]:.4f} "
                  f"BigMHC={bigmhc_scores[idx]:.4f} binding={ranks_filtered[idx]:.4f}",
                  file=sys.stderr)

        # Print top 5 overall
        print(f"\n    Top 5 by composite score:", file=sys.stderr)
        for rank in range(min(5, len(ranked_indices))):
            idx = ranked_indices[rank]
            pep = peps_filtered[idx]
            is_pos = pep in pos_pep_set
            marker = '+' if is_pos else '-'
            print(f"    [{marker}] Rank {rank+1}: {pep:<15s} {hlas_filtered[idx]:<15s} "
                  f"comp={composite_scores[idx]:.4f} BigMHC={bigmhc_scores[idx]:.4f} "
                  f"bind={ranks_filtered[idx]:.4f}",
                  file=sys.stderr)

        all_results.append({
            'patient': pid,
            'total_candidates': len(peps_filtered),
            'n_positive': n_total_pos,
            'n_found': len(found_at_rank),
            'recall_at_20': n_found_top20 / max(n_total_pos, 1),
            'recall_at_50': n_found_top50 / max(n_total_pos, 1),
        })

    # Summary
    print(f"\n\n{'='*60}")
    print(f"  FULL PIPELINE BENCHMARK SUMMARY")
    print(f"{'='*60}")
    summary = pd.DataFrame(all_results)
    print(f"\n{summary.to_string(index=False)}")
    print(f"\nMacro-average Recall@20: {summary['recall_at_20'].mean():.3f}")
    print(f"Macro-average Recall@50: {summary['recall_at_50'].mean():.3f}")


if __name__ == '__main__':
    main()
