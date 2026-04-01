#!/usr/bin/env python3
"""
Module 11: mRNA Design

Converts the polyepitope amino acid sequence into a synthesis-ready
mRNA nucleotide string:

  1. Codon optimization (LinearDesign if available, else CAI maximization)
  2. 5' UTR prepend
  3. 3' UTR append
  4. Poly-A tail
  5. Pseudouridine substitution flags
  6. Synthesis specification JSON

Inputs:
    polyepitope.faa         — from Module 10
    utr_5prime.fasta        — 5' UTR sequence
    utr_3prime.fasta        — 3' UTR sequence

Output:
    mrna_sequence.fasta     — final nucleotide string
    synthesis_spec.json     — specification for synthesis facility
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

# ── Human codon usage table (frequency per amino acid) ──
# Source: Kazusa codon usage database, Homo sapiens
HUMAN_CODON_TABLE = {
    'F': [('UUU', 0.45), ('UUC', 0.55)],
    'L': [('UUA', 0.07), ('UUG', 0.13), ('CUU', 0.13), ('CUC', 0.20), ('CUA', 0.07), ('CUG', 0.40)],
    'I': [('AUU', 0.36), ('AUC', 0.48), ('AUA', 0.16)],
    'M': [('AUG', 1.00)],
    'V': [('GUU', 0.18), ('GUC', 0.24), ('GUA', 0.11), ('GUG', 0.47)],
    'S': [('UCU', 0.15), ('UCC', 0.22), ('UCA', 0.12), ('UCG', 0.06), ('AGU', 0.15), ('AGC', 0.24)],  # removed duplicate
    'P': [('CCU', 0.28), ('CCC', 0.33), ('CCA', 0.27), ('CCG', 0.11)],
    'T': [('ACU', 0.24), ('ACC', 0.36), ('ACA', 0.28), ('ACG', 0.12)],
    'A': [('GCU', 0.26), ('GCC', 0.40), ('GCA', 0.23), ('GCG', 0.11)],
    'Y': [('UAU', 0.43), ('UAC', 0.57)],
    '*': [('UAA', 0.28), ('UAG', 0.20), ('UGA', 0.52)],
    'H': [('CAU', 0.41), ('CAC', 0.59)],
    'Q': [('CAA', 0.25), ('CAG', 0.75)],
    'N': [('AAU', 0.46), ('AAC', 0.54)],
    'K': [('AAA', 0.42), ('AAG', 0.58)],
    'D': [('GAU', 0.46), ('GAC', 0.54)],
    'E': [('GAA', 0.42), ('GAG', 0.58)],
    'C': [('UGU', 0.45), ('UGC', 0.55)],
    'W': [('UGG', 1.00)],
    'R': [('CGU', 0.08), ('CGC', 0.19), ('CGA', 0.11), ('CGG', 0.21), ('AGA', 0.20), ('AGG', 0.20)],
    'G': [('GGU', 0.16), ('GGC', 0.34), ('GGA', 0.25), ('GGG', 0.25)],
}


def read_fasta_sequence(fasta_path: str) -> str:
    """Read the first sequence from a FASTA file."""
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_parts:
                    break
                continue
            seq_parts.append(line)
    return ''.join(seq_parts)


def codon_optimize_cai(protein_seq: str) -> str:
    """
    Codon optimization using Codon Adaptation Index (CAI) maximization.
    Selects the most frequent human codon for each amino acid.
    """
    nt_parts = []
    for aa in protein_seq:
        if aa not in HUMAN_CODON_TABLE:
            raise ValueError(f"Unknown amino acid '{aa}' in protein sequence — cannot codon-optimize")
        codons = HUMAN_CODON_TABLE[aa]
        # Select highest-frequency codon
        best_codon = max(codons, key=lambda x: x[1])[0]
        nt_parts.append(best_codon)

    return ''.join(nt_parts)


def run_lineardesign(protein_seq: str, lineardesign_dir: str = None,
                      lam: float = 3.0) -> tuple[str, str, float, float]:
    """
    Run LinearDesign for joint codon + mRNA structure optimization.

    Calls the LinearDesign_2D binary directly (avoids Python 2 wrapper).
    On macOS ARM, runs via Docker (linux/amd64) since the pre-compiled
    .so is x86-only.

    Args:
        protein_seq: amino acid sequence
        lineardesign_dir: path to LinearDesign repo root
        lam: lambda parameter (0 = pure MFE, higher = more CAI weight)

    Returns:
        (mrna_sequence, structure, mfe_kcal, cai_score)

    Raises:
        RuntimeError if LinearDesign fails.
        FileNotFoundError if LinearDesign is not installed.
    """
    # Find LinearDesign directory
    search_paths = [
        lineardesign_dir,
        os.environ.get('LINEARDESIGN_DIR'),
        str(Path(__file__).parent.parent / 'LinearDesign'),
        '/opt/LinearDesign',
    ]

    ld_dir = None
    for p in search_paths:
        if p and Path(p).exists() and (Path(p) / 'bin' / 'LinearDesign_2D').exists():
            ld_dir = str(Path(p).resolve())
            break

    if ld_dir is None:
        # Check if we can run via Docker
        docker_check = subprocess.run(
            ['docker', 'info'], capture_output=True, timeout=5
        )
        if docker_check.returncode != 0:
            raise FileNotFoundError(
                "LinearDesign not found and Docker not available. "
                "Install LinearDesign or start Docker."
            )

        # Find the LinearDesign dir for Docker volume mount
        for p in search_paths:
            if p and Path(p).exists() and (Path(p) / 'Makefile').exists():
                ld_dir = str(Path(p).resolve())
                break

        if ld_dir is None:
            raise FileNotFoundError(
                "LinearDesign directory not found. Clone it: "
                "git clone https://github.com/LinearDesignSoftware/LinearDesign.git"
            )

    binary = Path(ld_dir) / 'bin' / 'LinearDesign_2D'
    codon_table = Path(ld_dir) / 'codon_usage_freq_table_human.csv'

    if binary.exists():
        # Try direct execution first
        try:
            result = subprocess.run(
                [str(binary), str(lam), '0', str(codon_table)],
                input=protein_seq,
                capture_output=True, text=True, timeout=600,
            )
            if result.returncode == 0:
                return _parse_lineardesign_output(result.stdout)
        except OSError:
            pass  # Binary incompatible — fall through to Docker

    # Run via Docker (needed on macOS ARM)
    print("  Running LinearDesign via Docker (linux/amd64)...", file=sys.stderr)
    result = subprocess.run(
        [
            'docker', 'run', '--rm', '--platform', 'linux/amd64',
            '-v', f'{ld_dir}:/build', '-w', '/build',
            '-i', 'gcc:11',
            'bash', '-c',
            f'echo "{protein_seq}" | ./bin/LinearDesign_2D {lam} 0 codon_usage_freq_table_human.csv',
        ],
        capture_output=True, text=True, timeout=600,
    )

    if result.returncode != 0:
        raise RuntimeError(f"LinearDesign failed: {result.stderr[:500]}")

    return _parse_lineardesign_output(result.stdout)


def _parse_lineardesign_output(stdout: str) -> tuple[str, str, float, float]:
    """Parse LinearDesign stdout into (mRNA_seq, structure, MFE, CAI)."""
    mrna_seq = None
    structure = None
    mfe = None
    cai = None

    for line in stdout.strip().split('\n'):
        line = line.strip()
        if line.startswith('mRNA sequence:'):
            mrna_seq = line.split(':', 1)[1].strip()
        elif line.startswith('mRNA structure:'):
            structure = line.split(':', 1)[1].strip()
        elif line.startswith('mRNA folding free energy:'):
            # "mRNA folding free energy: -35.70 kcal/mol; mRNA CAI: 0.898"
            parts = line.split(';')
            mfe_str = parts[0].split(':')[1].strip().split()[0]
            mfe = float(mfe_str)
            if len(parts) > 1:
                cai_str = parts[1].split(':')[1].strip()
                cai = float(cai_str)

    if mrna_seq is None:
        raise RuntimeError(f"Failed to parse LinearDesign output: {stdout[:300]}")

    # Validate mRNA sequence
    valid_nt = set('AUGC')
    if not all(c in valid_nt for c in mrna_seq):
        raise RuntimeError(f"LinearDesign produced invalid nucleotides: {mrna_seq[:50]}...")

    return mrna_seq, structure or '', mfe or 0.0, cai or 0.0


def count_pseudouridine_positions(mrna_seq: str) -> list[int]:
    """Find all uridine positions for N1-methylpseudouridine substitution."""
    return [i + 1 for i, nt in enumerate(mrna_seq) if nt == 'U']


def main():
    parser = argparse.ArgumentParser(description='Module 11: mRNA Design')
    parser.add_argument('--polyepitope', required=True, help='polyepitope.faa from Module 10')
    parser.add_argument('--utr-5prime', required=True, help='5\' UTR FASTA')
    parser.add_argument('--utr-3prime', required=True, help='3\' UTR FASTA')
    parser.add_argument('--poly-a-length', type=int, default=120)
    parser.add_argument('--patient-id', default='patient')
    parser.add_argument('--output-fasta', required=True)
    parser.add_argument('--lineardesign-dir', default=None,
                        help='Path to LinearDesign repo root')
    parser.add_argument('--lineardesign-lambda', type=float, default=3.0,
                        help='LinearDesign lambda (0=pure MFE, higher=more CAI)')
    parser.add_argument('--output-spec', required=True)
    args = parser.parse_args()

    # Read polyepitope protein sequence
    protein_seq = read_fasta_sequence(args.polyepitope)
    if not protein_seq:
        print("ERROR: Empty polyepitope sequence", file=sys.stderr)
        sys.exit(1)

    print(f"Input polyepitope: {len(protein_seq)} amino acids")

    # Read UTR sequences
    if not Path(args.utr_5prime).exists():
        raise FileNotFoundError(f"5' UTR file not found: {args.utr_5prime}")
    if not Path(args.utr_3prime).exists():
        raise FileNotFoundError(f"3' UTR file not found: {args.utr_3prime}")
    utr5 = read_fasta_sequence(args.utr_5prime)
    utr3 = read_fasta_sequence(args.utr_3prime)

    # ── Step 1: Codon optimization with LinearDesign ──
    mrna_result = run_lineardesign(
        protein_seq,
        lineardesign_dir=args.lineardesign_dir,
        lam=args.lineardesign_lambda,
    )
    cds, structure, mfe, cai = mrna_result
    codon_method = 'LinearDesign'
    print(f"Codon optimization: LinearDesign (lambda={args.lineardesign_lambda}, "
          f"MFE={mfe:.1f} kcal/mol, CAI={cai:.3f})")

    # Add stop codon
    stop_codon = 'UGA'  # most frequent human stop codon
    cds += stop_codon

    # ── Step 2: Assemble full mRNA ──
    # 5' cap is handled at synthesis (CleanCap-AG)
    mrna_parts = []
    if utr5:
        mrna_parts.append(utr5)
    mrna_parts.append(cds)
    if utr3:
        mrna_parts.append(utr3)

    poly_a = 'A' * args.poly_a_length
    mrna_parts.append(poly_a)

    mrna_sequence = ''.join(mrna_parts)

    # ── Step 3: Pseudouridine positions ──
    psi_positions = count_pseudouridine_positions(mrna_sequence)

    # ── Step 4: Write FASTA output ──
    with open(args.output_fasta, 'w') as f:
        f.write(f">{args.patient_id}_mrna len={len(mrna_sequence)} codon_opt={codon_method}\n")
        for i in range(0, len(mrna_sequence), 80):
            f.write(mrna_sequence[i:i + 80] + '\n')

    # ── Step 5: Write synthesis specification JSON ──
    # Compute neoantigen positions in nucleotide coordinates
    # CDS starts after 5' UTR
    cds_start = len(utr5) if utr5 else 0

    spec = {
        'patient_id': args.patient_id,
        'mrna_sequence': mrna_sequence,
        'length_nt': len(mrna_sequence),
        'cds_start': cds_start,
        'cds_end': cds_start + len(cds),
        'protein_length_aa': len(protein_seq),
        'codon_optimization': codon_method,
        'lineardesign_lambda': args.lineardesign_lambda,
        'lineardesign_mfe_kcal': mfe,
        'lineardesign_cai': cai,
        'lineardesign_structure': structure,
        'modifications': {
            'pseudouridine': True,
            'pseudouridine_type': 'N1-methylpseudouridine',
            'pseudouridine_count': len(psi_positions),
            'cap_analog': 'CleanCap-AG',
            'poly_a_length': args.poly_a_length,
        },
        'utr_5prime': utr5,
        'utr_5prime_length': len(utr5),
        'utr_3prime': utr3,
        'utr_3prime_length': len(utr3),
        'stop_codon': stop_codon,
        'gc_content': round(
            (mrna_sequence.count('G') + mrna_sequence.count('C')) / len(mrna_sequence), 4
        ),
        'synthesis_notes': (
            'N1-methylpseudouridine substitution at all uridine positions. '
            'CleanCap-AG co-transcriptional capping. '
            f'Poly-A tail: {args.poly_a_length}nt.'
        ),
    }

    with open(args.output_spec, 'w') as f:
        json.dump(spec, f, indent=2)

    # Summary
    print(f"\nmRNA design complete:")
    print(f"  Total length: {len(mrna_sequence)} nt")
    print(f"  5' UTR: {len(utr5)} nt")
    print(f"  CDS: {len(cds)} nt ({len(protein_seq)} aa + stop)")
    print(f"  3' UTR: {len(utr3)} nt")
    print(f"  Poly-A: {args.poly_a_length} nt")
    print(f"  GC content: {spec['gc_content']:.1%}")
    print(f"  Pseudouridine positions: {len(psi_positions)}")
    print(f"  Codon optimization: {codon_method}")


if __name__ == '__main__':
    main()
