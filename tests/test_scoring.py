"""Tests for Module 8a: Immunogenicity Scoring."""

import math
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'bin'))

from score_immunogenicity import (
    compute_agretopicity,
    compute_foreignness_kmer,
    build_proteome_kmers,
)


class TestAgretopicity:
    def test_frameshift_no_wildtype(self):
        score = compute_agretopicity(0.01, None)
        assert score == 5.0

    def test_mutation_creates_new_binder(self):
        score = compute_agretopicity(0.001, 0.5)
        assert score > 5.0

    def test_wildtype_already_binds(self):
        score = compute_agretopicity(0.01, 0.02)
        assert 0 < score < 2.0

    def test_wildtype_binds_better(self):
        score = compute_agretopicity(0.5, 0.01)
        assert score < 0

    def test_zero_wildtype(self):
        score = compute_agretopicity(0.01, 0)
        assert score == 5.0


class TestForeignness:
    def test_exact_match_returns_zero(self):
        kmers = {'YLQPRTFLL'}
        score = compute_foreignness_kmer('YLQPRTFLL', kmers, k=9)
        assert score == 0.0

    def test_no_match_returns_high(self):
        kmers = {'XXXXXXXXX'}
        score = compute_foreignness_kmer('YLQPRTFLL', kmers, k=9)
        assert score > 0.5

    def test_empty_proteome_raises(self):
        with pytest.raises(ValueError, match="Proteome k-mer set is empty"):
            compute_foreignness_kmer('YLQPRTFLL', set(), k=9)

    def test_short_peptide_raises(self):
        with pytest.raises(ValueError, match="shorter than k"):
            compute_foreignness_kmer('YLQ', {'YLQPRTFLL'}, k=9)


class TestProteomeKmers:
    def test_build_from_fasta(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">test_protein\n")
            f.write("MRGYLQPRTFLLKYFRGM\n")
            f.flush()
            kmers = build_proteome_kmers(f.name, k=9)
            assert 'YLQPRTFLL' in kmers
            assert len(kmers) == 10

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            build_proteome_kmers('/nonexistent/file.fasta', k=9)

    def test_none_path_raises(self):
        with pytest.raises(ValueError):
            build_proteome_kmers(None, k=9)
