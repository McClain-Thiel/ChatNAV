"""Tests for Module 11: mRNA Design."""

import json
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'bin'))

from design_mrna import (
    codon_optimize_cai,
    count_pseudouridine_positions,
    read_fasta_sequence,
    HUMAN_CODON_TABLE,
)


class TestCodonOptimization:
    def test_methionine(self):
        """M has only one codon: AUG."""
        result = codon_optimize_cai('M')
        assert result == 'AUG'

    def test_tryptophan(self):
        """W has only one codon: UGG."""
        result = codon_optimize_cai('W')
        assert result == 'UGG'

    def test_length(self):
        """Output should be exactly 3x the input protein length."""
        protein = 'METDTLLLWV'
        result = codon_optimize_cai(protein)
        assert len(result) == len(protein) * 3

    def test_valid_codons(self):
        """All output codons should be valid RNA (A, U, G, C)."""
        protein = 'ACDEFGHIKLMNPQRSTVWY'
        result = codon_optimize_cai(protein)
        assert all(nt in 'AUGC' for nt in result)

    def test_selects_most_frequent(self):
        """Should select the highest-frequency codon for each AA."""
        for aa, codons in HUMAN_CODON_TABLE.items():
            if aa == '*':
                continue
            best = max(codons, key=lambda x: x[1])[0]
            result = codon_optimize_cai(aa)
            assert result == best, f"For AA '{aa}', expected {best}, got {result}"


class TestPseudouridine:
    def test_count_positions(self):
        seq = 'AUGCUUAGC'
        positions = count_pseudouridine_positions(seq)
        # U at positions 2, 5, 6 (1-indexed): A-U-G-C-U-U-A-G-C
        assert 2 in positions
        assert 5 in positions
        assert 6 in positions
        assert len(positions) == 3

    def test_no_uridines(self):
        seq = 'AGGCCCAAA'
        positions = count_pseudouridine_positions(seq)
        assert len(positions) == 0

    def test_all_uridines(self):
        seq = 'UUUUU'
        positions = count_pseudouridine_positions(seq)
        assert positions == [1, 2, 3, 4, 5]


class TestReadFasta:
    def test_single_sequence(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">test\n")
            f.write("METDTLLL\n")
            f.write("WVLLLWVP\n")
            f.flush()

            seq = read_fasta_sequence(f.name)
            assert seq == 'METDTLLLWVLLLWVP'

    def test_multiline(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write("AAAA\n")
            f.write("BBBB\n")
            f.write(">seq2\n")
            f.write("CCCC\n")
            f.flush()

            seq = read_fasta_sequence(f.name)
            assert seq == 'AAAABBBB'  # only first sequence


class TestSynthesisSpec:
    def test_full_pipeline(self):
        """Integration test: protein → codon opt → UTR → poly-A → spec."""
        protein = 'METDTLLLWVLLLWVPGSTGDYLQPRTFLL'

        # Codon optimize
        cds = codon_optimize_cai(protein) + 'UGA'

        # UTRs
        utr5 = 'GAAAUAAGAGAGAAAAGAAGAG'
        utr3 = 'UGAUAAUAGGCUGGAGCCUCGG'

        # Assemble
        mrna = utr5 + cds + utr3 + 'A' * 120

        assert len(mrna) > 0
        assert mrna.startswith(utr5)
        assert mrna.endswith('A' * 120)
        assert len(cds) == len(protein) * 3 + 3  # +3 for stop codon

        # GC content should be reasonable (40-70%)
        gc = (mrna.count('G') + mrna.count('C')) / len(mrna)
        assert 0.3 < gc < 0.8

        # Pseudouridine count should be > 0
        psi = count_pseudouridine_positions(mrna)
        assert len(psi) > 0
