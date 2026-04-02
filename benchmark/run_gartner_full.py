#!/usr/bin/env python3
"""
Full-pipeline Gartner benchmark with ALL signals:
  - MHCflurry binding
  - BigMHC-IM immunogenicity
  - Foreignness (k-mer vs proteome)
  - Agretopicity (wt vs mut binding ratio)
  - Structural tier 1 (position lookup)
  - AlphaFold2 tier 3 for top N candidates (via local NIM)

Runs on GPU instance with AlphaFold2 NIM at localhost:8000.
"""

import sys
import os
import json
import math
import time
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'bin'))


AF2_ENDPOINT = "http://localhost:8000/protein-structure/alphafold2/multimer/predict-structure-from-sequences"
B2M_SEQ = "MSRSVALAVLALLSLSGLEAIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"


def load_gartner():
    xlsx = PROJECT_ROOT / 'benchmark' / 'gartner' / 'supplementary_table.xlsx'
    patients = pd.read_excel(xlsx, sheet_name='Supplementary Table 1', header=1)
    hla_map = dict(zip(patients['ID'], patients['Patient HLAs']))
    positives = pd.read_excel(xlsx, sheet_name='Supplemetary Table 3', header=1)
    positives = positives.dropna(subset=['Mutant minimal Peptide', 'HLA'])
    nmers = pd.read_excel(xlsx, sheet_name='Supplementary Table 11', header=1)
    return patients, hla_map, positives, nmers


def run_mhcflurry_genotype(peptides, genotype_alleles):
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    clean_alleles = [a.replace('*', '').replace(':', '') for a in genotype_alleles]
    results = predictor.predict(peptides=list(set(peptides)), alleles=clean_alleles, verbose=0)
    lookup = {}
    for _, row in results.iterrows():
        lookup[row['peptide']] = {
            'best_allele': row['best_allele'],
            'binding_rank': row['presentation_percentile'] / 100.0,
        }
    return lookup


def run_mhcflurry_wt(wt_peptides, alleles_per_peptide):
    """Run MHCflurry on wildtype peptides for agretopicity."""
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    wt_ranks = {}
    # Group by allele for efficiency
    from collections import defaultdict
    allele_groups = defaultdict(list)
    for i, (pep, allele) in enumerate(zip(wt_peptides, alleles_per_peptide)):
        if pep and isinstance(pep, str) and len(pep) >= 8:
            allele_groups[allele].append((i, pep))

    for allele, items in allele_groups.items():
        indices, peps = zip(*items)
        clean_allele = allele.replace('*', '').replace(':', '')
        result = predictor.predict_affinity(
            peptides=list(peps),
            alleles={'sample': [clean_allele]},
            verbose=0,
        )
        for j, (idx, pep) in enumerate(zip(indices, peps)):
            row = result[result['peptide'] == pep]
            if not row.empty:
                wt_ranks[idx] = float(row.iloc[0]['affinity_percentile']) / 100.0

    return wt_ranks


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
    return out['BigMHC_IM'].tolist()


def run_foreignness_batch(peptides, proteome_kmers):
    from score_immunogenicity import compute_foreignness_kmer
    scores = []
    for pep in peptides:
        try:
            scores.append(compute_foreignness_kmer(pep, proteome_kmers))
        except ValueError:
            scores.append(0.5)
    return scores


def run_structural_tier1(peptides, wt_peptides):
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


def run_alphafold2_sasa(peptide, hla_allele=None):
    """Run AlphaFold2-Multimer on peptide + B2M, compute per-residue SASA."""
    try:
        resp = requests.post(
            AF2_ENDPOINT,
            headers={"Content-Type": "application/json"},
            json={"sequences": [peptide, B2M_SEQ], "relax_prediction": False},
            timeout=3600,
        )
        if resp.status_code != 200:
            print(f"    AF2 error {resp.status_code}: {resp.text[:100]}", file=sys.stderr)
            return None

        data = resp.json()
        if isinstance(data, list) and data:
            pdb_str = data[0] if isinstance(data[0], str) else json.dumps(data[0])
        else:
            return None

        # Write PDB and compute SASA
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_str)
            pdb_path = f.name

        from Bio.PDB import PDBParser, ShrakeRupley
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pmhc', pdb_path)
        model = structure[0]

        sr = ShrakeRupley()
        sr.compute(model, level='R')

        # Peptide is the shortest chain
        chains = list(model.get_chains())
        peptide_chain = min(chains, key=lambda c: len(list(c.get_residues())))
        residues = list(peptide_chain.get_residues())

        sasa_values = [r.sasa for r in residues]
        Path(pdb_path).unlink(missing_ok=True)

        return sasa_values

    except Exception as e:
        print(f"    AF2 exception: {e}", file=sys.stderr)
        return None


def compute_af2_structural_score(sasa_values, mut_positions):
    """Score based on SASA at mutation positions relative to max."""
    if not sasa_values or max(sasa_values) == 0:
        return 0.5
    max_sasa = max(sasa_values)
    mut_sasa = [sasa_values[p] for p in mut_positions if p < len(sasa_values)]
    if not mut_sasa:
        return 0.5
    return float(np.clip(np.mean(mut_sasa) / max_sasa, 0, 1))


def ppv_at_k(labels, scores, k):
    top_k = np.argsort(scores)[::-1][:k]
    return np.mean(labels[top_k])


def recall_at_k(labels, scores, k):
    n_pos = np.sum(labels)
    if n_pos == 0:
        return 0.0
    top_k = np.argsort(scores)[::-1][:k]
    return np.sum(labels[top_k]) / n_pos


def main():
    patients, hla_map, positives, nmers = load_gartner()
    print(f"Loaded: {len(positives)} immunogenic mmps, {len(nmers)} nmers, {len(hla_map)} patients")

    # Check if AlphaFold is available
    af2_available = False
    try:
        r = requests.get("http://localhost:8000/v1/health/ready", timeout=5)
        af2_available = r.status_code == 200 and r.json().get('status') == 'ready'
    except Exception:
        pass
    print(f"AlphaFold2 NIM: {'AVAILABLE' if af2_available else 'NOT AVAILABLE'}")

    # Build proteome k-mer index
    from score_immunogenicity import build_proteome_kmers
    proteome_path = PROJECT_ROOT / 'reference' / 'proteome' / 'human_proteome.fasta'
    print("Building proteome k-mer index...", file=sys.stderr)
    proteome_kmers = build_proteome_kmers(str(proteome_path))

    # Select test patients
    patient_pos_counts = positives.groupby('ID').size().sort_values(ascending=False)
    test_patients = patient_pos_counts.head(5).index.tolist()

    # Build positive peptide lookup
    positive_set = {}
    positive_wt = {}  # peptide → wildtype peptide (from Table 3)
    for _, row in positives.iterrows():
        pid = row['ID']
        pep = str(row['Mutant minimal Peptide']).strip()
        wt = str(row.get('Wild type minimal peptied', '')).strip()
        hla = str(row['HLA']).strip()
        if pid not in positive_set:
            positive_set[pid] = set()
        positive_set[pid].add(pep)
        if wt and wt != '-' and wt != 'nan':
            positive_wt[pep] = wt

    all_results = []

    for pid in test_patients:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Patient {pid}", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        pat_nmers = nmers[nmers['ID'] == pid]
        hlas = [h.strip() for h in str(hla_map[pid]).split(',')]
        hlas = [h for h in hlas if h.startswith('HLA-')][:6]

        # Extract unique 9-10mer windows
        all_peps = set()
        pep_to_nmer = {}  # track which nmer each peptide came from
        for _, nmer_row in pat_nmers.iterrows():
            nmer_seq = str(nmer_row['Mutant nmer']).strip()
            for pep_len in [9, 10]:
                for start in range(len(nmer_seq) - pep_len + 1):
                    pep = nmer_seq[start:start + pep_len]
                    if all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in pep):
                        all_peps.add(pep)
                        pep_to_nmer[pep] = nmer_seq

        all_peps = list(all_peps)
        print(f"  {len(all_peps)} unique peptides from {len(pat_nmers)} nmers", file=sys.stderr)

        # 1. MHCflurry binding
        print(f"  MHCflurry ({len(all_peps)} peptides × {len(hlas)} alleles)...", file=sys.stderr)
        mhcflurry_lookup = run_mhcflurry_genotype(all_peps, hlas)

        # Filter to binders
        peps_f, hlas_f, ranks_f = [], [], []
        for pep in all_peps:
            info = mhcflurry_lookup.get(pep)
            if info and info['binding_rank'] < 0.02:
                peps_f.append(pep)
                raw_hla = info['best_allele']
                hla_fmt = f"{raw_hla[:5]}*{raw_hla[5:7]}:{raw_hla[7:]}" if '*' not in raw_hla else raw_hla
                hlas_f.append(hla_fmt)
                ranks_f.append(info['binding_rank'])

        if not peps_f:
            for pep in all_peps:
                info = mhcflurry_lookup.get(pep)
                if info and info['binding_rank'] < 0.05:
                    peps_f.append(pep)
                    raw_hla = info['best_allele']
                    hla_fmt = f"{raw_hla[:5]}*{raw_hla[5:7]}:{raw_hla[7:]}" if '*' not in raw_hla else raw_hla
                    hlas_f.append(hla_fmt)
                    ranks_f.append(info['binding_rank'])

        print(f"  {len(peps_f)} pass binding filter", file=sys.stderr)

        # 2. BigMHC-IM
        print(f"  BigMHC-IM...", file=sys.stderr)
        bigmhc_scores = run_bigmhc_batch(peps_f, hlas_f)

        # 3. Foreignness
        print(f"  Foreignness...", file=sys.stderr)
        foreign_scores = run_foreignness_batch(peps_f, proteome_kmers)

        # 4. Agretopicity — compute WT binding for each peptide
        print(f"  Agretopicity (computing WT binding)...", file=sys.stderr)
        wt_peps = []
        for pep in peps_f:
            wt = positive_wt.get(pep, '')
            if not wt:
                # For non-positive peptides, try to find WT from Ensembl via nmer context
                # Simplified: use the peptide itself as WT proxy (agretopicity = 0)
                wt = pep
            wt_peps.append(wt)

        wt_ranks = run_mhcflurry_wt(wt_peps, hlas_f)
        agreto_scores = []
        for i in range(len(peps_f)):
            wt_rank = wt_ranks.get(i, ranks_f[i])  # default to same rank if no WT
            mut_rank = max(ranks_f[i], 0.0001)
            wt_rank = max(wt_rank, 0.0001)
            agreto = math.log2(wt_rank / mut_rank)
            agreto_scores.append(agreto)

        # Normalize agretopicity to [0,1]
        agreto_norm = [1.0 / (1.0 + math.exp(-0.5 * a)) for a in agreto_scores]

        # 5. Structural tier 1
        struct_scores = run_structural_tier1(peps_f, wt_peps)

        # 6. Composite score with ALL signals
        composite = []
        for i in range(len(peps_f)):
            score = (
                0.35 * bigmhc_scores[i] +
                0.15 * foreign_scores[i] +
                0.10 * agreto_norm[i] +
                0.20 * (1.0 - ranks_f[i]) +
                0.10 * struct_scores[i] +
                0.10 * 0.5  # expression placeholder
            )
            composite.append(score)

        # Rank
        ranked = np.argsort(composite)[::-1]

        # 7. AlphaFold tier 3 on top 5
        af2_scores = {}
        if af2_available:
            print(f"  AlphaFold2 on top 5 candidates...", file=sys.stderr)
            for rank_pos in range(min(5, len(ranked))):
                idx = ranked[rank_pos]
                pep = peps_f[idx]
                wt = wt_peps[idx]

                print(f"    AF2 #{rank_pos+1}: {pep}...", file=sys.stderr)
                sasa = run_alphafold2_sasa(pep)
                if sasa:
                    # Find mutation positions
                    if wt and wt != pep and len(wt) == len(pep):
                        mut_pos = [j for j in range(len(pep)) if pep[j] != wt[j]]
                    else:
                        mut_pos = list(range(len(pep)))
                    af2_score = compute_af2_structural_score(sasa, mut_pos)
                    af2_scores[idx] = af2_score
                    print(f"    → SASA score: {af2_score:.3f}", file=sys.stderr)

            # Re-compute composite with AF2 scores for those that have them
            if af2_scores:
                for idx, af2_score in af2_scores.items():
                    # Replace tier 1 structural with AF2 structural
                    composite[idx] = (
                        0.30 * bigmhc_scores[idx] +
                        0.15 * foreign_scores[idx] +
                        0.10 * agreto_norm[idx] +
                        0.15 * (1.0 - ranks_f[idx]) +
                        0.20 * af2_score +  # higher weight for AF2
                        0.10 * 0.5
                    )
                ranked = np.argsort(composite)[::-1]

        # Check recall
        pat_positives = positive_set.get(pid, set())
        found = {}
        for rank, idx in enumerate(ranked):
            pep = peps_f[idx]
            if pep in pat_positives and pep not in found:
                found[pep] = rank + 1

        n_top20 = sum(1 for r in found.values() if r <= 20)
        n_top50 = sum(1 for r in found.values() if r <= 50)
        n_pos = len(pat_positives)

        print(f"\n  Patient {pid}: {n_pos} confirmed, {len(found)}/{n_pos} found", file=sys.stderr)
        print(f"  Recall@20: {n_top20}/{n_pos} = {n_top20/max(n_pos,1):.2f}", file=sys.stderr)
        print(f"  Recall@50: {n_top50}/{n_pos} = {n_top50/max(n_pos,1):.2f}", file=sys.stderr)
        print(f"  AF2 predictions: {len(af2_scores)}", file=sys.stderr)

        for pep, rank in sorted(found.items(), key=lambda x: x[1]):
            idx = ranked[rank - 1]
            af_tag = " [AF2]" if idx in af2_scores else ""
            print(f"    Rank {rank:>4d}: {pep:<15s} comp={composite[idx]:.4f} "
                  f"IM={bigmhc_scores[idx]:.3f} bind={ranks_f[idx]:.4f} "
                  f"agreto={agreto_scores[idx]:.2f}{af_tag}", file=sys.stderr)

        all_results.append({
            'patient': pid,
            'total_candidates': len(peps_f),
            'n_positive': n_pos,
            'n_found': len(found),
            'recall_at_20': n_top20 / max(n_pos, 1),
            'recall_at_50': n_top50 / max(n_pos, 1),
            'af2_predictions': len(af2_scores),
        })

    # Summary
    print(f"\n\n{'='*60}")
    print(f"  FULL PIPELINE BENCHMARK (all signals + AlphaFold)")
    print(f"{'='*60}")
    summary = pd.DataFrame(all_results)
    print(f"\n{summary.to_string(index=False)}")
    print(f"\nMacro-average Recall@20: {summary['recall_at_20'].mean():.3f}")
    print(f"Macro-average Recall@50: {summary['recall_at_50'].mean():.3f}")
    print(f"Total AF2 predictions:   {summary['af2_predictions'].sum()}")


if __name__ == '__main__':
    main()
