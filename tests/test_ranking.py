"""Tests for Module 9: Candidate Ranking and Selection."""

import math
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / 'bin'))

from rank_and_select import (
    normalize_expression,
    load_filter_profile,
    apply_sequential_filters,
    apply_rank_bonuses,
    ensure_hla_diversity,
)

WEIGHTS_YAML = str(Path(__file__).parent.parent / 'conf' / 'scoring_weights.yaml')


class TestNormalization:
    def test_expression_zero(self):
        assert normalize_expression(0) == pytest.approx(0.0, abs=0.01)

    def test_expression_high(self):
        # TPM of 10000 → log10(10001)/5 ≈ 0.8
        val = normalize_expression(10000)
        assert 0.7 < val < 0.9

    def test_expression_cap(self):
        # Should not exceed 1.0
        val = normalize_expression(1e10)
        assert val <= 1.0

    def test_expression_negative(self):
        assert normalize_expression(-5) == 0.0

    def test_expression_none(self):
        assert normalize_expression(None) == 0.0


class TestFilterProfile:
    def test_load_high_tmb(self):
        profile = load_filter_profile(WEIGHTS_YAML, 'high_tmb')
        assert 'hard_filters' in profile
        assert 'bonuses' in profile
        assert profile['hard_filters']['mhc_binding_rank'] == 0.02

    def test_load_low_tmb(self):
        profile = load_filter_profile(WEIGHTS_YAML, 'low_tmb')
        assert profile['bonuses']['is_frameshift'] == 0.7

    def test_load_research(self):
        profile = load_filter_profile(WEIGHTS_YAML, 'research')
        assert profile['hard_filters']['mhc_binding_rank'] == 1.0

    def test_unknown_profile_raises(self):
        with pytest.raises(ValueError, match="not found"):
            load_filter_profile(WEIGHTS_YAML, 'nonexistent')


class TestSequentialFilters:
    def _sample_df(self):
        return pd.DataFrame({
            'peptide_id': ['p1', 'p2', 'p3', 'p4'],
            'binding_rank': [0.001, 0.05, 0.01, 0.001],
            'tpm': [50.0, 0.1, 5.0, 100.0],
            'ccf': [0.9, 0.8, 0.3, 0.7],
        })

    def test_binding_rank_filter(self):
        df = self._sample_df()
        filters = {'mhc_binding_rank': 0.02, 'min_tpm': 0, 'min_ccf': 0}
        result = apply_sequential_filters(df, filters)
        assert len(result) == 3  # p2 filtered (rank 0.05)

    def test_tpm_filter(self):
        df = self._sample_df()
        filters = {'mhc_binding_rank': 1.0, 'min_tpm': 1.0, 'min_ccf': 0}
        result = apply_sequential_filters(df, filters)
        assert len(result) == 3  # p2 filtered (tpm 0.1)

    def test_ccf_filter(self):
        df = self._sample_df()
        filters = {'mhc_binding_rank': 1.0, 'min_tpm': 0, 'min_ccf': 0.5}
        result = apply_sequential_filters(df, filters)
        assert len(result) == 3  # p3 filtered (ccf 0.3)

    def test_combined_filters(self):
        df = self._sample_df()
        filters = {'mhc_binding_rank': 0.02, 'min_tpm': 1.0, 'min_ccf': 0.5}
        result = apply_sequential_filters(df, filters)
        assert len(result) == 2  # p1 and p4 pass


class TestRankBonuses:
    def test_frameshift_promotion(self):
        df = pd.DataFrame({
            'binding_rank': [0.01, 0.01],
            'tpm': [10.0, 10.0],
            'is_frameshift': [True, False],
        })
        bonuses = {'is_frameshift': 0.5}
        result = apply_rank_bonuses(df, bonuses)
        # Frameshift should have lower effective_rank (= better)
        assert result.iloc[0]['effective_rank'] < result.iloc[1]['effective_rank']


class TestHLADiversity:
    def test_already_diverse(self):
        df = pd.DataFrame({
            'composite_score': [0.9, 0.8, 0.7],
            'hla_allele': ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-A*02:01'],
        })
        result = ensure_hla_diversity(df, 2)
        assert len(result) == 2

    def test_forces_diversity(self):
        df = pd.DataFrame({
            'composite_score': [0.9, 0.85, 0.8, 0.7],
            'hla_allele': ['HLA-A*02:01', 'HLA-A*02:01', 'HLA-A*02:01', 'HLA-B*07:02'],
        })
        result = ensure_hla_diversity(df, 3)
        assert len(result) == 3
        assert result['hla_allele'].nunique() >= 2
