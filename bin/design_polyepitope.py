#!/usr/bin/env python3
"""
Module 10: Polyepitope Design

Assembles selected neoantigens into an optimized polyepitope construct:
  1. Separate MHC-I and MHC-II epitopes
  2. Optimize ordering to minimize junctional epitope formation (greedy TSP)
  3. Insert appropriate linkers (AAY for MHC-I, GPGPG for MHC-II)
  4. Prepend signal peptide for secretory pathway / MHC-II loading
  5. Optionally append MITD domain
  6. QC: check junction regions for unintended strong MHC binders

Inputs:
    selected_neoantigens.tsv    — from Module 9
    signal_peptides.fasta       — signal peptide + MITD sequences

Output:
    polyepitope.faa             — full polyepitope amino acid sequence
    polyepitope_design.tsv      — epitope order, linkers, junctional QC
"""

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

# ── Constants ──
LINKER_MHC_I = 'AAY'
LINKER_MHC_II = 'GPGPG'
TERMINAL_LINKER = 'AAAAA'
SIGNAL_PEPTIDE = 'METDTLLLWVLLLWVPGSTGD'

# Amino acids with high MHC-I anchor potential (simplified)
ANCHOR_RESIDUES = set('LFIVMYW')


def load_signal_sequences(fasta_path: str) -> dict:
    """Load signal peptide and MITD from FASTA file."""
    seqs = {}
    current_name = None
    current_seq = []

    if not Path(fasta_path).exists():
        raise FileNotFoundError(f"Signal peptide file not found: {fasta_path}")

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    seqs[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0].lower()
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            seqs[current_name] = ''.join(current_seq)

    return seqs


def junction_score(pep_a: str, pep_b: str, linker: str, k: int = 9) -> float:
    """
    Score a junction region for unintended MHC-I epitope formation.

    Extracts all k-mers spanning the junction (last few AAs of pep_a +
    linker + first few AAs of pep_b) and scores them based on anchor
    residue patterns. Lower score = safer junction.

    This is a simplified heuristic; in production, call NetMHCpan on
    junction peptides.
    """
    junction_seq = pep_a[-4:] + linker + pep_b[:4]
    score = 0.0

    for i in range(len(junction_seq) - k + 1):
        kmer = junction_seq[i:i + k]
        # Simple heuristic: anchor residues at positions 2 and 9 (C-terminus)
        # indicate potential MHC-I binding
        if len(kmer) >= 9:
            pos2 = kmer[1] in ANCHOR_RESIDUES
            pos9 = kmer[-1] in ANCHOR_RESIDUES
            if pos2 and pos9:
                score += 1.0
            elif pos2 or pos9:
                score += 0.3

    return score


def greedy_tsp_ordering(peptides: list[str], linker: str) -> list[int]:
    """
    Greedy TSP heuristic: order peptides to minimize total junction score.

    Start with the peptide that has the lowest average junction score
    with all others. At each step, append the unvisited peptide with
    the lowest junction score to the current tail.
    """
    n = len(peptides)
    if n <= 2:
        return list(range(n))

    # Precompute pairwise junction scores
    scores = {}
    for i, j in itertools.permutations(range(n), 2):
        scores[(i, j)] = junction_score(peptides[i], peptides[j], linker)

    # Start with the peptide that has lowest average outgoing score
    avg_scores = []
    for i in range(n):
        avg = sum(scores.get((i, j), 0) for j in range(n) if j != i) / (n - 1)
        avg_scores.append((avg, i))
    avg_scores.sort()

    order = [avg_scores[0][1]]
    remaining = set(range(n)) - {order[0]}

    while remaining:
        current = order[-1]
        best_next = min(remaining, key=lambda j: scores.get((current, j), float('inf')))
        order.append(best_next)
        remaining.remove(best_next)

    return order


def build_polyepitope(epitopes_i: list[dict], epitopes_ii: list[dict],
                       signal_seq: str, mitd_seq: str) -> tuple[str, list[dict]]:
    """
    Build the full polyepitope construct:
        signal_peptide + [MHC-I epitopes with AAY linkers] +
        [MHC-II epitopes with GPGPG linkers] + terminal_linker + MITD
    """
    design_rows = []
    construct_parts = []

    # Signal peptide
    construct_parts.append(signal_seq)
    design_rows.append({
        'position': 0,
        'element': 'signal_peptide',
        'sequence': signal_seq,
        'linker': None,
        'junction_score': None,
    })

    # MHC-I epitopes
    if epitopes_i:
        peptides_i = [e['peptide_sequence'] for e in epitopes_i]
        order_i = greedy_tsp_ordering(peptides_i, LINKER_MHC_I)

        for idx, pos in enumerate(order_i):
            ep = epitopes_i[pos]
            pep = ep['peptide_sequence']

            if idx > 0:
                # Compute junction score with previous epitope
                prev_pep = epitopes_i[order_i[idx - 1]]['peptide_sequence']
                jscore = junction_score(prev_pep, pep, LINKER_MHC_I)
                construct_parts.append(LINKER_MHC_I)
            else:
                jscore = None

            construct_parts.append(pep)
            design_rows.append({
                'position': idx + 1,
                'element': f"mhc_i_epitope",
                'peptide_id': ep.get('peptide_id', ''),
                'gene': ep.get('gene', ''),
                'mutation': ep.get('mutation', ''),
                'rank': ep.get('rank', ''),
                'sequence': pep,
                'hla_allele': ep.get('hla_allele', ''),
                'linker': LINKER_MHC_I if idx > 0 else None,
                'junction_score': round(jscore, 3) if jscore is not None else None,
            })

    # Transition linker between MHC-I and MHC-II sections
    if epitopes_i and epitopes_ii:
        construct_parts.append(LINKER_MHC_II)

    # MHC-II epitopes
    if epitopes_ii:
        peptides_ii = [e['peptide_sequence'] for e in epitopes_ii]
        order_ii = greedy_tsp_ordering(peptides_ii, LINKER_MHC_II)

        for idx, pos in enumerate(order_ii):
            ep = epitopes_ii[pos]
            pep = ep['peptide_sequence']

            if idx > 0:
                prev_pep = epitopes_ii[order_ii[idx - 1]]['peptide_sequence']
                jscore = junction_score(prev_pep, pep, LINKER_MHC_II)
                construct_parts.append(LINKER_MHC_II)
            else:
                jscore = None

            construct_parts.append(pep)
            design_rows.append({
                'position': len(epitopes_i) + idx + 1 if epitopes_i else idx + 1,
                'element': f"mhc_ii_epitope",
                'peptide_id': ep.get('peptide_id', ''),
                'gene': ep.get('gene', ''),
                'mutation': ep.get('mutation', ''),
                'rank': ep.get('rank', ''),
                'sequence': pep,
                'hla_allele': ep.get('hla_allele', ''),
                'linker': LINKER_MHC_II if idx > 0 else None,
                'junction_score': round(jscore, 3) if jscore is not None else None,
            })

    # Terminal linker
    construct_parts.append(TERMINAL_LINKER)

    # MITD domain (optional)
    if mitd_seq:
        construct_parts.append(mitd_seq)
        design_rows.append({
            'position': len(design_rows),
            'element': 'mitd',
            'sequence': mitd_seq,
            'linker': TERMINAL_LINKER,
            'junction_score': None,
        })

    polyepitope_seq = ''.join(construct_parts)
    return polyepitope_seq, design_rows


def main():
    parser = argparse.ArgumentParser(description='Module 10: Polyepitope Design')
    parser.add_argument('--selected-neoantigens', required=True)
    parser.add_argument('--signal-peptides', required=True)
    parser.add_argument('--output-fasta', required=True)
    parser.add_argument('--output-design', required=True)
    parser.add_argument('--patient-id', default='patient')
    args = parser.parse_args()

    # Load selected neoantigens
    df = pd.read_csv(args.selected_neoantigens, sep='\t')

    if len(df) == 0:
        print("ERROR: No selected neoantigens to design polyepitope from", file=sys.stderr)
        sys.exit(1)

    # Load signal sequences
    signal_seqs = load_signal_sequences(args.signal_peptides)
    signal_peptide = signal_seqs.get('signal_peptide', SIGNAL_PEPTIDE)
    mitd = signal_seqs.get('mitd', '')

    # Separate MHC-I and MHC-II epitopes
    epitopes_i = []
    epitopes_ii = []

    for _, row in df.iterrows():
        ep = row.to_dict()
        # Require peptide_sequence — no fallbacks
        if 'peptide_sequence' not in ep or pd.isna(ep.get('peptide_sequence')):
            raise ValueError(
                f"Peptide {ep.get('peptide_id', '?')} has no peptide_sequence. "
                f"selected_neoantigens.tsv must include peptide_sequence column."
            )

        mhc_class = ep.get('mhc_class')
        if mhc_class is None or pd.isna(mhc_class):
            raise ValueError(f"Peptide {ep.get('peptide_id', '?')} has no mhc_class — required for polyepitope design")
        mhc_class = str(mhc_class)
        if mhc_class == 'II':
            epitopes_ii.append(ep)
        else:
            epitopes_i.append(ep)

    print(f"Designing polyepitope: {len(epitopes_i)} MHC-I + {len(epitopes_ii)} MHC-II epitopes")

    # Build the construct
    polyepitope_seq, design_rows = build_polyepitope(
        epitopes_i, epitopes_ii, signal_peptide, mitd
    )

    # Write FASTA output
    with open(args.output_fasta, 'w') as f:
        f.write(f">{args.patient_id}_polyepitope len={len(polyepitope_seq)}\n")
        # Write in 80-char lines
        for i in range(0, len(polyepitope_seq), 80):
            f.write(polyepitope_seq[i:i + 80] + '\n')

    # Write design TSV
    design_df = pd.DataFrame(design_rows)
    design_df.to_csv(args.output_design, sep='\t', index=False)

    # Summary
    print(f"Polyepitope: {len(polyepitope_seq)} amino acids")
    print(f"  Signal peptide: {len(signal_peptide)} aa")
    print(f"  MHC-I epitopes: {len(epitopes_i)}")
    print(f"  MHC-II epitopes: {len(epitopes_ii)}")
    if mitd:
        print(f"  MITD domain: {len(mitd)} aa")

    # Flag high junction scores
    high_junctions = [r for r in design_rows
                      if r.get('junction_score') is not None and r['junction_score'] > 1.0]
    if high_junctions:
        print(f"\nWARNING: {len(high_junctions)} junctions with score > 1.0 — review for junctional epitopes",
              file=sys.stderr)
        for j in high_junctions:
            print(f"  Position {j['position']}: {j.get('gene', '?')} — score {j['junction_score']}",
                  file=sys.stderr)


if __name__ == '__main__':
    main()
