#!/usr/bin/env python3
"""
Module 8a: Immunogenicity Scoring

Computes three immunogenicity signals for each (peptide, HLA) pair:

  1. BigMHC-IM immunogenicity score (Nature Machine Intelligence 2023)
  2. Foreignness score — 1 - max(similarity to human proteome)
  3. Agretopicity — log2(wt_binding_rank / mut_binding_rank)

Inputs:
  binding_predictions.tsv   — from Module 7 / Step 0
  candidates_meta.tsv       — from Module 6 / Step 0
  human_proteome.fasta      — UniProt reference proteome

Output:
  immunogenicity_scores.tsv
"""

import argparse
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

def _run_bigmhc_inprocess(input_path: str):
    """Run BigMHC-IM prediction in the current process."""
    import torch
    from bigmhc.src.cli import parseArgs
    from bigmhc.src.bigmhc import BigMHC

    original_argv = sys.argv
    try:
        sys.argv = [
            'bigmhc_predict',
            '-i', input_path,
            '-m', 'im',
            '-p', '0',
            '-a', '1',
            '-c', '0',
            '-d', 'cpu',
        ]
        args, data, models = parseArgs(train=False)

        # Predict (from BigMHC predict.py entry point)
        for x in range(len(models)):
            models[x].eval()
            models[x] = BigMHC.accelerate(models[x], devices=args.devices)

        preds = []
        with torch.no_grad():
            for idx, bat in enumerate(data):
                out = []
                for model in models:
                    dev = next(model.parameters()).device
                    _out, _att = model(
                        mhc=bat.mhc.to(dev),
                        pep=bat.pep.to(dev))
                    out.append(torch.sigmoid(_out))
                out = torch.mean(torch.stack(out), dim=0)
                rawbat = data.dataset.getbat(idx=idx, enc=False)
                rawbat[args.modelname] = out.cpu().numpy()
                preds.append(rawbat)

        result_df = pd.concat(preds).sort_index()
        result_df.to_csv(input_path + '.prd', index=False)
    finally:
        sys.argv = original_argv


def predict_bigmhc_batch(peptides: list[str], hla_alleles: list[str]) -> list[float]:
    """
    Predict immunogenicity using BigMHC-IM.

    BigMHC-IM is a transformer-based model trained on immunogenicity data.
    Supports 8-14mer peptides. Returns scores in [0, 1] where higher = more immunogenic.

    Raises RuntimeError if BigMHC fails.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        # Write peptide,HLA pairs
        for pep, hla in zip(peptides, hla_alleles):
            f.write(f"{pep},{hla}\n")
        input_path = f.name

    output_path = input_path + '.prd'

    try:
        # Run BigMHC in-process
        _run_bigmhc_inprocess(input_path)

        if not Path(output_path).exists():
            raise RuntimeError(f"BigMHC produced no output file at {output_path}")

        # Parse output
        out_df = pd.read_csv(output_path)
        if 'BigMHC_IM' not in out_df.columns:
            raise RuntimeError(f"BigMHC output missing BigMHC_IM column. Columns: {list(out_df.columns)}")

        scores = out_df['BigMHC_IM'].tolist()

        if len(scores) != len(peptides):
            raise RuntimeError(
                f"BigMHC returned {len(scores)} scores for {len(peptides)} peptides"
            )

        return [float(s) for s in scores]

    finally:
        Path(input_path).unlink(missing_ok=True)
        Path(output_path).unlink(missing_ok=True)


def compute_foreignness_kmer(peptide: str, proteome_kmers: set, k: int = 9) -> float:
    """
    Foreignness score: how dissimilar is this peptide to anything in the
    human self-proteome?

    Uses a k-mer overlap approach: for each k-mer in the peptide, check
    how many 1-mismatch matches exist in the proteome k-mer set.

    Returns a score in [0, 1] where 1 = maximally foreign (no close match).
    """
    if not proteome_kmers:
        raise ValueError("Proteome k-mer set is empty — provide --human-proteome")
    if len(peptide) < k:
        raise ValueError(f"Peptide '{peptide}' is shorter than k={k} — cannot compute foreignness")

    best_similarity = 0.0
    for i in range(len(peptide) - k + 1):
        kmer = peptide[i:i + k]
        if kmer in proteome_kmers:
            return 0.0  # exact match in self → not foreign at all

        # Count 1-mismatch matches
        for pos in range(k):
            for aa in 'ACDEFGHIKLMNPQRSTVWY':
                variant = kmer[:pos] + aa + kmer[pos + 1:]
                if variant in proteome_kmers:
                    similarity = (k - 1) / k
                    best_similarity = max(best_similarity, similarity)

    return 1.0 - best_similarity


def build_proteome_kmers(proteome_path: str, k: int = 9) -> set:
    """Build k-mer set from a FASTA proteome file."""
    if not proteome_path:
        raise ValueError("proteome_path is required — provide a human proteome FASTA")
    if not Path(proteome_path).exists():
        raise FileNotFoundError(f"Proteome file not found: {proteome_path}")

    kmers = set()
    seq = []
    with open(proteome_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    protein = ''.join(seq)
                    for i in range(len(protein) - k + 1):
                        kmers.add(protein[i:i + k])
                seq = []
            else:
                seq.append(line)
        if seq:
            protein = ''.join(seq)
            for i in range(len(protein) - k + 1):
                kmers.add(protein[i:i + k])

    return kmers


def compute_agretopicity(mut_rank: float, wt_rank: float) -> float:
    """
    Agretopicity: how much better does the mutant bind vs wildtype?

    agretopicity = log2(wt_rank / mut_rank)

    High positive value = mutation creates a new binding event.
    Zero or negative = wildtype already binds well → tolerance risk.

    For frameshifts (no wildtype), returns a high default value.
    """
    if wt_rank is None or pd.isna(wt_rank) or wt_rank <= 0:
        return 5.0  # frameshift / no wildtype → maximally foreign

    if mut_rank <= 0:
        mut_rank = 0.0001  # avoid log(0)

    return math.log2(wt_rank / mut_rank)


def main():
    parser = argparse.ArgumentParser(description='Module 8a: Immunogenicity Scoring')
    parser.add_argument('--binding-predictions', required=True,
                        help='binding_predictions.tsv from MHCflurry')
    parser.add_argument('--candidates-meta', required=True,
                        help='candidates_meta.tsv with peptide sequences')
    parser.add_argument('--human-proteome', required=True,
                        help='Human reference proteome FASTA for foreignness scoring')
    parser.add_argument('--output', required=True,
                        help='Output immunogenicity_scores.tsv')
    args = parser.parse_args()

    # Load inputs
    binding = pd.read_csv(args.binding_predictions, sep='\t')
    meta = pd.read_csv(args.candidates_meta, sep='\t')

    # Build proteome k-mer index for foreignness scoring
    proteome_kmers = build_proteome_kmers(args.human_proteome)
    print(f"Built {len(proteome_kmers)} k-mers from human proteome", file=sys.stderr)

    # Merge binding predictions with candidate metadata
    merged = binding.merge(meta, on='peptide_id', how='left')

    # Validate required columns
    if 'peptide_sequence' not in merged.columns:
        raise ValueError("candidates_meta.tsv must contain 'peptide_sequence' column")
    missing_seqs = merged['peptide_sequence'].isna().sum()
    if missing_seqs > 0:
        raise ValueError(f"{missing_seqs} peptides have no peptide_sequence")

    peptide_seqs = merged['peptide_sequence'].tolist()
    hla_alleles = merged['hla_allele'].tolist()

    # Batch predict immunogenicity with BigMHC-IM
    print(f"Running BigMHC-IM on {len(peptide_seqs)} peptide-HLA pairs...", file=sys.stderr)
    immuno_scores = predict_bigmhc_batch(peptide_seqs, hla_alleles)
    print(f"BigMHC-IM scoring complete", file=sys.stderr)

    results = []
    for idx, (_, row) in enumerate(merged.iterrows()):
        peptide = row['peptide_id']
        peptide_seq = peptide_seqs[idx]
        hla = row['hla_allele']

        mut_rank = row.get('binding_rank')
        if mut_rank is None or pd.isna(mut_rank):
            raise ValueError(f"Peptide {peptide} has no binding_rank — MHC binding data required")

        is_frameshift = row.get('is_frameshift', False)

        # 1. BigMHC-IM immunogenicity score
        immuno = immuno_scores[idx]

        # 2. Foreignness (skip for peptides shorter than k=9)
        if len(peptide_seq) >= 9:
            foreignness = compute_foreignness_kmer(peptide_seq, proteome_kmers)
        else:
            foreignness = 0.5  # 8-mers can't be k-mer matched at k=9

        # 3. Agretopicity
        wt_rank = row.get('wildtype_binding_rank', None)
        if (wt_rank is None or (isinstance(wt_rank, float) and pd.isna(wt_rank)) or wt_rank == '') and not is_frameshift:
            raise ValueError(
                f"Peptide {peptide} is not a frameshift but has no wildtype_binding_rank. "
                f"Provide wildtype binding data or mark as frameshift."
            )

        agretopicity = compute_agretopicity(mut_rank, wt_rank)

        results.append({
            'peptide_id': peptide,
            'hla_allele': hla,
            'immunogenicity_score': round(immuno, 5),
            'foreignness_score': round(foreignness, 5),
            'agretopicity': round(agretopicity, 5),
        })

    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(out_df)} immunogenicity scores to {args.output}")


if __name__ == '__main__':
    main()
