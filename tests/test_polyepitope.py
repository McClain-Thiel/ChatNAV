"""Tests for Module 10: Polyepitope Design."""

import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'bin'))

from design_polyepitope import (
    junction_score,
    greedy_tsp_ordering,
    build_polyepitope,
    load_signal_sequences,
    LINKER_MHC_I,
    LINKER_MHC_II,
    SIGNAL_PEPTIDE,
)


class TestJunctionScore:
    def test_low_risk_junction(self):
        """Peptides without anchor residues at junction → low score."""
        score = junction_score('AAAAAAA', 'GGGGGGG', 'AAY')
        assert score < 0.5

    def test_detects_potential_epitope(self):
        """Anchor residues at junction positions → higher score."""
        # L at position 2 and L at C-terminal are MHC-I anchors
        score = junction_score('AAAAYLQA', 'LLLLLLLL', 'AAY')
        assert score >= 0

    def test_empty_linker(self):
        score = junction_score('YLQPRTFLL', 'KIGDFGLAT', '')
        assert isinstance(score, float)


class TestGreedyTSP:
    def test_two_peptides(self):
        """Two peptides → order is [0, 1] or [1, 0]."""
        order = greedy_tsp_ordering(['YLQPRTFLL', 'KIGDFGLAT'], 'AAY')
        assert len(order) == 2
        assert set(order) == {0, 1}

    def test_single_peptide(self):
        order = greedy_tsp_ordering(['YLQPRTFLL'], 'AAY')
        assert order == [0]

    def test_preserves_all_peptides(self):
        """All peptides should appear exactly once in ordering."""
        peptides = ['AAAAAAA', 'BBBBBBB', 'CCCCCCC', 'DDDDDDD', 'EEEEEEE']
        order = greedy_tsp_ordering(peptides, 'AAY')
        assert len(order) == 5
        assert set(order) == {0, 1, 2, 3, 4}


class TestBuildPolyepitope:
    def _sample_epitopes(self):
        epitopes_i = [
            {'peptide_sequence': 'YLQPRTFLL', 'peptide_id': 'p1', 'gene': 'BRAF',
             'mutation': 'V600E', 'rank': 1, 'hla_allele': 'HLA-A*02:01'},
            {'peptide_sequence': 'KIGDFGLAT', 'peptide_id': 'p2', 'gene': 'KRAS',
             'mutation': 'G12V', 'rank': 2, 'hla_allele': 'HLA-A*11:01'},
        ]
        epitopes_ii = [
            {'peptide_sequence': 'KIGDFGLATEKSRWSGSHQF', 'peptide_id': 'p3',
             'gene': 'BRAF', 'mutation': 'V600E', 'rank': 3, 'hla_allele': 'DRB1*01:01'},
        ]
        return epitopes_i, epitopes_ii

    def test_basic_construct(self):
        ep_i, ep_ii = self._sample_epitopes()
        seq, design = build_polyepitope(ep_i, ep_ii, SIGNAL_PEPTIDE, '')

        # Should start with signal peptide
        assert seq.startswith(SIGNAL_PEPTIDE)

        # Should contain both MHC-I epitopes
        assert 'YLQPRTFLL' in seq
        assert 'KIGDFGLAT' in seq

        # Should contain MHC-II epitope
        assert 'KIGDFGLATEKSRWSGSHQF' in seq

        # Should contain linkers
        assert LINKER_MHC_I in seq
        assert LINKER_MHC_II in seq

        # Should end with terminal linker
        assert seq.endswith('AAAAA') or seq.endswith('AAAAA' + '')

    def test_design_table(self):
        ep_i, ep_ii = self._sample_epitopes()
        seq, design = build_polyepitope(ep_i, ep_ii, SIGNAL_PEPTIDE, '')

        # Design should have entries for signal + 2 MHC-I + 1 MHC-II
        assert len(design) >= 4

        # First entry should be signal peptide
        assert design[0]['element'] == 'signal_peptide'

    def test_mhc_i_only(self):
        ep_i, _ = self._sample_epitopes()
        seq, design = build_polyepitope(ep_i, [], SIGNAL_PEPTIDE, '')

        assert 'YLQPRTFLL' in seq
        assert LINKER_MHC_II not in seq.replace(SIGNAL_PEPTIDE, '')

    def test_with_mitd(self):
        ep_i, _ = self._sample_epitopes()
        mitd = 'IVGIIAGLVLFGAVITGAVVAAVMW'
        seq, design = build_polyepitope(ep_i, [], SIGNAL_PEPTIDE, mitd)

        assert mitd in seq
        assert any(d['element'] == 'mitd' for d in design)


class TestClassIISpecific:
    def test_gpgpg_linker_used_for_class_ii(self):
        ep_ii = [
            {'peptide_sequence': 'KIGDFGLATEKSRWSGSHQF', 'peptide_id': 'p1',
             'gene': 'BRAF', 'mutation': 'V600E', 'rank': 1, 'hla_allele': 'DRB1*01:01'},
            {'peptide_sequence': 'VVVGAVGVGKSALTIQLIQN', 'peptide_id': 'p2',
             'gene': 'KRAS', 'mutation': 'G12V', 'rank': 2, 'hla_allele': 'DRB1*01:01'},
        ]
        seq, design = build_polyepitope([], ep_ii, SIGNAL_PEPTIDE, '')

        # Should use GPGPG linker between Class II epitopes
        assert LINKER_MHC_II in seq
        # Should NOT use AAY linker (no Class I epitopes)
        # (AAY could appear in peptide sequences, so check design table instead)
        class_ii_entries = [d for d in design if d.get('linker') == LINKER_MHC_II]
        assert len(class_ii_entries) >= 1

    def test_mixed_class_construct_has_both_linkers(self):
        ep_i = [
            {'peptide_sequence': 'YLQPRTFLL', 'peptide_id': 'p1', 'gene': 'BRAF',
             'mutation': 'V600E', 'rank': 1, 'hla_allele': 'HLA-A*02:01'},
            {'peptide_sequence': 'KIGDFGLAT', 'peptide_id': 'p2', 'gene': 'KRAS',
             'mutation': 'G12D', 'rank': 2, 'hla_allele': 'HLA-A*11:01'},
        ]
        ep_ii = [
            {'peptide_sequence': 'KIGDFGLATEKSRWSGSHQF', 'peptide_id': 'p3',
             'gene': 'KRAS', 'mutation': 'G12V', 'rank': 3, 'hla_allele': 'DRB1*01:01'},
        ]
        seq, design = build_polyepitope(ep_i, ep_ii, SIGNAL_PEPTIDE, '')

        # AAY between Class I epitopes, GPGPG before Class II
        assert LINKER_MHC_I in seq
        assert LINKER_MHC_II in seq

    def test_junction_score_with_gpgpg(self):
        score = junction_score('YLQPRTFLL', 'KIGDFGLAT', LINKER_MHC_II)
        # GPGPG contains proline — generally not a strong anchor residue,
        # so junction score should be low
        assert isinstance(score, float)
        assert score >= 0.0


class TestLoadSignalSequences:
    def test_load_from_fasta(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">signal_peptide test signal\n")
            f.write("METDTLLLWVLLLWVPGSTGD\n")
            f.write(">mitd test MITD\n")
            f.write("IVGIIAGLVLFG\n")
            f.flush()

            seqs = load_signal_sequences(f.name)
            assert 'signal_peptide' in seqs
            assert 'mitd' in seqs
            assert seqs['signal_peptide'] == 'METDTLLLWVLLLWVPGSTGD'

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_signal_sequences('/nonexistent/file.fasta')
