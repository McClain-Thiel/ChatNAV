#!/usr/bin/env python3
"""
Module 8b: Structural Scoring (Tiered)

Three-tier approach to score whether the mutated residue is TCR-exposed:

  Tier 1 — Position lookup (instant, free)
    For canonical 9-10mer MHC-I peptides with known alleles:
    Is the mutation at a TCR-facing position (P3-P7 for 9-mers)?
    Returns a binary flag + confidence score based on position.

  Tier 2 — PANDORA local modeling (1-2 min, CPU, free)
    For long peptides (>10mer) or unusual alleles where position
    rules are unreliable. Runs PANDORA homology modeling to generate
    a pMHC structure, then computes per-residue solvent-accessible
    surface area (SASA) at the mutation site.
    Requires: MODELLER license (free academic), PANDORA installed.

  Tier 3 — AlphaFold3 via NVIDIA API (expensive, top 5-10 only)
    Full pMHC structure prediction for the highest-scoring candidates.
    Reserved for final validation. Requires: NVIDIA API key + GPU server.

Inputs:
    immunogenicity_scores.tsv
    candidates_meta.tsv (with peptide_sequence and wildtype_peptide)

Output:
    structural_scores.tsv — structural_score + tier used for each candidate
"""

import argparse
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd
import requests

# ── Tier 1: Position lookup tables ──

# TCR-facing positions for MHC-I peptides (0-indexed)
# Based on crystal structure surveys (Trolle et al. 2016, Rasmussen et al. 2016)
# These positions point "up" toward the TCR in the MHC groove
TCR_FACING_POSITIONS = {
    8:  {2, 3, 4, 5},           # 8-mer: P3-P6
    9:  {3, 4, 5, 6},           # 9-mer: P4-P7 (canonical)
    10: {3, 4, 5, 6, 7},        # 10-mer: P4-P8
    11: {3, 4, 5, 6, 7, 8},     # 11-mer: P4-P9
}

# Primary anchor positions (buried in MHC groove, NOT TCR-facing)
ANCHOR_POSITIONS = {
    8:  {1, 7},     # P2, P-omega
    9:  {1, 8},     # P2, P9
    10: {1, 9},     # P2, P10
    11: {1, 10},    # P2, P11
}

# Position weights — central positions are more exposed to TCR
# Based on average SASA from crystal structure surveys
POSITION_WEIGHTS_9MER = {
    0: 0.1,   # P1 — partially buried
    1: 0.0,   # P2 — anchor (buried)
    2: 0.3,   # P3 — partially exposed
    3: 0.7,   # P4 — highly exposed
    4: 0.9,   # P5 — most exposed (peak of bulge)
    5: 0.8,   # P6 — highly exposed
    6: 0.5,   # P7 — partially exposed
    7: 0.2,   # P8 — partially buried
    8: 0.0,   # P9 — anchor (buried)
}


def find_mutation_positions(mut_peptide: str, wt_peptide: str) -> list[int]:
    """Find 0-indexed positions where mutant differs from wildtype."""
    positions = []
    for i, (m, w) in enumerate(zip(mut_peptide, wt_peptide)):
        if m != w:
            positions.append(i)
    return positions


def score_tier1_position(peptide: str, mut_positions: list[int],
                          peptide_length: int) -> tuple[float, str]:
    """
    Tier 1: Score based on known TCR-facing positions.

    Returns (score, rationale) where score is in [0, 1].
    Higher = mutation more likely TCR-exposed.
    """
    if peptide_length not in TCR_FACING_POSITIONS:
        return None, f"No position rules for {peptide_length}-mer"

    tcr_positions = TCR_FACING_POSITIONS[peptide_length]
    anchor_positions = ANCHOR_POSITIONS[peptide_length]

    if not mut_positions:
        return 0.0, "No mutation detected in peptide"

    # Check each mutation position
    scores = []
    details = []
    for pos in mut_positions:
        if pos in anchor_positions:
            scores.append(0.0)
            details.append(f"P{pos+1}=anchor(buried)")
        elif pos in tcr_positions:
            # Use position weight for 9-mers, interpolate for others
            if peptide_length == 9 and pos in POSITION_WEIGHTS_9MER:
                weight = POSITION_WEIGHTS_9MER[pos]
            else:
                # For non-9-mers, central positions get higher weight
                center = peptide_length / 2
                dist_from_center = abs(pos - center) / center
                weight = max(0, 1.0 - dist_from_center)
            scores.append(weight)
            details.append(f"P{pos+1}=TCR-facing({weight:.1f})")
        else:
            scores.append(0.15)
            details.append(f"P{pos+1}=edge({0.15})")

    avg_score = np.mean(scores)
    rationale = "; ".join(details)
    return float(avg_score), rationale


def score_tier2_pandora(peptide: str, hla_allele: str, mut_positions: list[int],
                         output_dir: str) -> tuple[float, str]:
    """
    Tier 2: PANDORA homology modeling + SASA computation.

    Builds a pMHC structure with PANDORA, then computes solvent-accessible
    surface area at the mutation position(s) using BioPython ShrakeRupley.

    Returns (score, rationale) where score is relative SASA at mutation site.
    """
    try:
        from PANDORA import Target, Pandora, Database
        from Bio.PDB import ShrakeRupley
    except ImportError:
        raise ImportError(
            "PANDORA not installed. Install with: pip install CSB-PANDORA\n"
            "Also requires MODELLER license: https://salilab.org/modeller/registration.html"
        )

    # Normalize HLA format for PANDORA: HLA-A*02:01 → HLA-A*0201
    hla_pandora = hla_allele.replace(':', '')

    # Load template database
    db = Database.load()

    # Define target
    pep_len = len(peptide)
    anchors = [2, pep_len]  # P2 and P-omega (canonical MHC-I anchors)

    target = Target(
        id=f'structural_{peptide}',
        allele_type=hla_pandora,
        peptide=peptide,
        anchors=anchors,
        MHC_class='I',
        output_dir=output_dir,
    )

    # Run modeling (fast mode: 1 model, very fast refinement)
    case = Pandora.Pandora(target, db)
    case.model(n_loop_models=1, loop_refinement='very_fast')

    if not case.results:
        raise RuntimeError(f"PANDORA produced no models for {peptide}/{hla_allele}")

    # Get best model
    best = case.results[0]
    structure = best.pdb
    model = structure[0]

    # Compute SASA
    sr = ShrakeRupley()
    sr.compute(model, level='R')

    # Extract SASA for peptide chain (chain P in PANDORA output)
    peptide_chain = model['P']
    residues = list(peptide_chain.get_residues())

    all_sasa = [r.sasa for r in residues]
    mut_sasa = [residues[pos].sasa for pos in mut_positions if pos < len(residues)]

    if not all_sasa or max(all_sasa) == 0:
        return 0.5, "SASA computation failed (all zeros)"

    # Normalize: mutation SASA relative to max SASA in peptide
    max_sasa = max(all_sasa)
    avg_mut_sasa = np.mean(mut_sasa) if mut_sasa else 0
    score = avg_mut_sasa / max_sasa

    sasa_str = ", ".join(f"P{p+1}={residues[p].sasa:.1f}" for p in mut_positions if p < len(residues))
    rationale = f"PANDORA SASA: {sasa_str}; max={max_sasa:.1f}; relative={score:.3f}"

    return float(np.clip(score, 0, 1)), rationale


def submit_tier3_alphafold(peptide_ids: list[str], hla_seqs: dict,
                            peptides: list[str], hla_alleles: list[str],
                            api_key: str) -> str:
    """
    Tier 3: Submit AlphaFold batch to Tamarind.bio API.

    Submits pMHC complexes (HLA heavy chain + B2M + peptide) as a batch.
    Returns batch ID for polling.
    """
    B2M = (
        "MSRSVALAVLALLSLSGLEAIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVE"
        "HSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
    )

    base_url = "https://app.tamarind.bio/api/"
    headers = {'x-api-key': api_key}

    settings = []
    job_names = []
    for pep_id, peptide, hla in zip(peptide_ids, peptides, hla_alleles):
        # Look up HLA heavy chain sequence
        hla_key = hla.replace(':', '').replace('*', '*')
        hla_seq = hla_seqs.get(hla, hla_seqs.get(hla_key, ''))
        if not hla_seq:
            raise ValueError(f"No HLA heavy chain sequence for {hla}")

        # AlphaFold multimer: concatenate chains with ':'
        multimer_seq = f"{hla_seq}:{B2M}:{peptide}"
        settings.append({"sequence": multimer_seq})
        job_names.append(pep_id)

    params = {
        "batchName": "neoantigen_structural_scoring",
        "type": "alphafold",
        "settings": settings,
        "jobNames": job_names,
    }

    resp = requests.post(
        base_url + "submit-batch",
        headers=headers,
        json=params,
        timeout=60,
    )

    if resp.status_code != 200:
        raise RuntimeError(f"Tamarind batch submit failed: {resp.status_code} {resp.text[:300]}")

    result = resp.json()
    batch_id = result.get('batchId') or result.get('id')
    print(f"  Submitted AlphaFold batch: {batch_id} ({len(settings)} jobs)", file=sys.stderr)
    return batch_id


def poll_tier3_results(batch_id: str, api_key: str,
                        timeout: int = 7200, interval: int = 60) -> list[dict]:
    """
    Poll Tamarind.bio for batch results.
    Returns list of completed job results.
    """
    base_url = "https://app.tamarind.bio/api/"
    headers = {'x-api-key': api_key}
    elapsed = 0

    while elapsed < timeout:
        resp = requests.get(
            base_url + f"batch/{batch_id}",
            headers=headers,
            timeout=30,
        )

        if resp.status_code != 200:
            raise RuntimeError(f"Batch status check failed: {resp.status_code}")

        data = resp.json()
        status = data.get('status', '')

        completed = sum(1 for j in data.get('jobs', []) if j.get('status') == 'completed')
        total = len(data.get('jobs', []))
        print(f"  Batch {batch_id}: {completed}/{total} complete ({elapsed}s elapsed)",
              file=sys.stderr)

        if status in ('completed', 'done', 'finished'):
            return data.get('jobs', [])

        if all(j.get('status') in ('completed', 'failed') for j in data.get('jobs', [])):
            return data.get('jobs', [])

        time.sleep(interval)
        elapsed += interval

    raise TimeoutError(f"AlphaFold batch {batch_id} did not complete within {timeout}s")


def score_tier3_from_pdb(pdb_string: str, mut_peptide: str, wt_peptide: str) -> tuple[float, str]:
    """
    Score TCR exposure from AlphaFold PDB output using SASA.
    The peptide is the last chain in the multimer prediction.
    """
    from Bio.PDB import PDBParser, ShrakeRupley
    import io

    mut_positions = find_mutation_positions(mut_peptide, wt_peptide) if wt_peptide else list(range(len(mut_peptide)))

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pmhc', io.StringIO(pdb_string))
    model = structure[0]

    # Peptide is the last/shortest chain
    chains = list(model.get_chains())
    peptide_chain = min(chains, key=lambda c: len(list(c.get_residues())))
    residues = list(peptide_chain.get_residues())

    sr = ShrakeRupley()
    sr.compute(model, level='R')

    all_sasa = [r.sasa for r in residues]
    mut_sasa = [residues[pos].sasa for pos in mut_positions if pos < len(residues)]

    if not all_sasa or max(all_sasa) == 0:
        return 0.5, "SASA all zeros from AF3 PDB"

    max_sasa = max(all_sasa)
    avg_mut_sasa = np.mean(mut_sasa) if mut_sasa else 0
    score = avg_mut_sasa / max_sasa

    sasa_str = ", ".join(f"P{p+1}={residues[p].sasa:.1f}" for p in mut_positions if p < len(residues))
    return float(np.clip(score, 0, 1)), f"AF3 SASA: {sasa_str}; max={max_sasa:.1f}; rel={score:.3f}"


def main():
    parser = argparse.ArgumentParser(description='Module 8b: Structural Scoring (Tiered)')
    parser.add_argument('--immunogenicity-scores', required=True)
    parser.add_argument('--candidates-meta', required=True)
    parser.add_argument('--enable-tier2', action='store_true',
                        help='Enable PANDORA modeling for ambiguous cases (requires MODELLER)')
    parser.add_argument('--enable-tier3', action='store_true',
                        help='Enable AlphaFold3 via Tamarind.bio for top candidates')
    parser.add_argument('--tier3-top-n', type=int, default=10,
                        help='Number of top candidates to send to AlphaFold3')
    parser.add_argument('--hla-heavy-chains', default=None,
                        help='HLA heavy chain FASTA (required for tier 3)')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    scores_df = pd.read_csv(args.immunogenicity_scores, sep='\t')
    meta_df = pd.read_csv(args.candidates_meta, sep='\t')

    # Merge to get peptide sequences and wildtype
    merged = scores_df.merge(
        meta_df[['peptide_id', 'peptide_sequence', 'wildtype_peptide']],
        on='peptide_id', how='left',
    )

    results = []
    tier_counts = {'tier1': 0, 'tier2': 0, 'tier3': 0, 'skipped': 0}

    with tempfile.TemporaryDirectory(prefix='structural_') as tmpdir:
        for _, row in merged.iterrows():
            peptide_id = row['peptide_id']
            hla = row['hla_allele']
            mut_pep = row.get('peptide_sequence', '')
            wt_pep = row.get('wildtype_peptide', '')

            if not mut_pep or pd.isna(mut_pep):
                results.append({
                    'peptide_id': peptide_id,
                    'hla_allele': hla,
                    'structural_score': None,
                    'structural_tier': 'skipped',
                    'structural_rationale': 'No peptide sequence',
                })
                tier_counts['skipped'] += 1
                continue

            # Find mutation positions
            if wt_pep and not pd.isna(wt_pep) and len(wt_pep) == len(mut_pep):
                mut_positions = find_mutation_positions(mut_pep, wt_pep)
            else:
                mut_positions = list(range(len(mut_pep)))

            pep_len = len(mut_pep)

            # Tier 1: Position lookup
            score, rationale = score_tier1_position(mut_pep, mut_positions, pep_len)

            if score is not None:
                tier = 'tier1'
                tier_counts['tier1'] += 1
            elif args.enable_tier2:
                try:
                    score, rationale = score_tier2_pandora(
                        mut_pep, hla, mut_positions, tmpdir
                    )
                    tier = 'tier2'
                    tier_counts['tier2'] += 1
                except Exception as e:
                    score = 0.5
                    rationale = f"Tier 2 failed: {e}"
                    tier = 'tier2_failed'
                    tier_counts['skipped'] += 1
            else:
                score = 0.5
                rationale = f"No position rules for {pep_len}-mer; tier2 not enabled"
                tier = 'tier1_unsupported'
                tier_counts['skipped'] += 1

            results.append({
                'peptide_id': peptide_id,
                'hla_allele': hla,
                'structural_score': round(score, 5) if score is not None else None,
                'structural_tier': tier,
                'structural_rationale': rationale,
            })

    # ── Tier 3: AlphaFold via Tamarind.bio for top candidates ──
    if args.enable_tier3:
        tamarind_key = os.environ.get('TAMARIND_API_KEY')
        if not tamarind_key:
            raise RuntimeError(
                "TAMARIND_API_KEY environment variable required for tier 3. "
                "Get a key at https://app.tamarind.bio/"
            )
        if not args.hla_heavy_chains:
            raise RuntimeError("--hla-heavy-chains required for tier 3 AlphaFold")

        # Load HLA sequences
        hla_seqs = {}
        with open(args.hla_heavy_chains) as f:
            name, parts = None, []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if name:
                        hla_seqs[name] = ''.join(parts)
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line)
            if name:
                hla_seqs[name] = ''.join(parts)

        # Sort by tier1 score, take top N
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('structural_score', ascending=False)
        top_n = results_df.head(args.tier3_top_n)

        # Get peptide sequences for top candidates
        meta_lookup = meta_df.set_index('peptide_id')
        tier3_ids, tier3_peps, tier3_hlas, tier3_wts = [], [], [], []
        for _, row in top_n.iterrows():
            pid = row['peptide_id']
            if pid in meta_lookup.index:
                meta_row = meta_lookup.loc[pid]
                pep_seq = meta_row['peptide_sequence'] if isinstance(meta_row, pd.Series) else meta_row.iloc[0]['peptide_sequence']
                wt_seq = meta_row.get('wildtype_peptide', '') if isinstance(meta_row, pd.Series) else ''
                tier3_ids.append(pid)
                tier3_peps.append(pep_seq)
                tier3_hlas.append(row['hla_allele'])
                tier3_wts.append(wt_seq if isinstance(wt_seq, str) else '')

        print(f"\nTier 3: Submitting top {len(tier3_ids)} candidates to AlphaFold...",
              file=sys.stderr)

        batch_id = submit_tier3_alphafold(
            tier3_ids, hla_seqs, tier3_peps, tier3_hlas, tamarind_key
        )

        # Poll for results
        jobs = poll_tier3_results(batch_id, tamarind_key)

        # Score completed jobs
        for job in jobs:
            if job.get('status') != 'completed':
                continue
            job_name = job.get('jobName', '')
            pdb_data = job.get('result', {}).get('pdb', '') or job.get('pdbString', '')
            if not pdb_data or not job_name:
                continue

            # Find matching result and update
            idx = next((i for i, r in enumerate(results) if r['peptide_id'] == job_name), None)
            if idx is not None:
                wt = tier3_wts[tier3_ids.index(job_name)] if job_name in tier3_ids else ''
                pep = tier3_peps[tier3_ids.index(job_name)]
                score, rationale = score_tier3_from_pdb(pdb_data, pep, wt)
                results[idx]['structural_score'] = round(score, 5)
                results[idx]['structural_tier'] = 'tier3'
                results[idx]['structural_rationale'] = rationale
                tier_counts['tier3'] += 1

    # Output
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep='\t', index=False)

    print(f"\nStructural scoring: {len(results)} candidates", file=sys.stderr)
    print(f"  Tier 1 (position lookup): {tier_counts['tier1']}", file=sys.stderr)
    print(f"  Tier 2 (PANDORA):         {tier_counts['tier2']}", file=sys.stderr)
    print(f"  Tier 3 (AlphaFold):       {tier_counts['tier3']}", file=sys.stderr)
    print(f"  Skipped:                  {tier_counts['skipped']}", file=sys.stderr)
    print(f"  → {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
