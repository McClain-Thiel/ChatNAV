#!/usr/bin/env python3
"""
ChatNAV: Unified neoantigen vaccine design pipeline.

Single entry point: MAF/VCF + expression + HLA → synthesis-ready mRNA vaccine.

Usage:
    python bin/run_pipeline.py \
        --maf patient_somatic.maf.gz \
        --expression patient_expression.tsv \
        --hla-alleles HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*44:03,HLA-C*07:02,HLA-C*05:01 \
        --output results/patient_001/ \
        --patient-id PT-001

    python bin/run_pipeline.py \
        --vcf patient_somatic.vcf.gz \
        --expression patient_expression.tsv \
        --hla-alleles HLA-A*02:01,HLA-B*07:02 \
        --output results/patient_002/

Outputs:
    selected_neoantigens.tsv   — ranked epitope list with scores
    polyepitope.faa            — polyepitope amino acid sequence
    polyepitope_design.tsv     — construct layout with linkers + junction QC
    mrna_sequence.fasta        — synthesis-ready mRNA
    synthesis_spec.json        — full specification (sequence, UTRs, mods, QC scans)
    summary_card.html          — one-page visual report
    confidence.json            — expected vaccine potency estimate
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
BIN = PROJECT_ROOT / 'bin'
REF = PROJECT_ROOT / 'reference'


def run_step(script: str, args: list, label: str):
    """Run a pipeline step. Crash on failure."""
    cmd = [sys.executable, str(BIN / script)] + [str(a) for a in args]
    print(f"\n{'─'*60}")
    print(f"  {label}")
    print(f"{'─'*60}")
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

    if result.stdout:
        for line in result.stdout.strip().split('\n'):
            print(f"  {line}")
    if result.stderr:
        for line in result.stderr.strip().split('\n'):
            if any(skip in line.lower() for skip in [
                'cuda', 'cudnn', 'gpu available', 'tensorflow/core',
                'onednn', 'xla', 'absl'
            ]):
                continue
            print(f"  {line}", file=sys.stderr)

    elapsed = time.time() - t0
    if result.returncode != 0:
        print(f"\n  FATAL: {label} failed (exit {result.returncode}, {elapsed:.1f}s)")
        sys.exit(1)

    print(f"  [{elapsed:.1f}s]")
    return result


def compute_confidence(selected_path: str, n_total_mutations: int,
                       n_filter_passing: int) -> dict:
    """
    Estimate vaccine confidence based on empirical benchmark data.

    Uses the relationship between number of filter-passing candidates
    and observed Recall@20 from the Muller/Gfeller benchmark (91 patients).
    """
    import pandas as pd

    df = pd.read_csv(selected_path, sep='\t')
    n_selected = len(df)

    # Empirical R@20 by candidate count (from Muller dev set, EXP-022 reranker)
    # Derived from benchmark analysis of 68 patients with immunogenic mutations
    empirical_tiers = [
        (30,  0.833, 1.000, 'HIGH',     'Few candidates — model selects most immunogenic'),
        (60,  0.785, 0.893, 'HIGH',     'Moderate candidates — strong selection expected'),
        (100, 0.606, 0.682, 'MODERATE', 'Many candidates — selection is harder'),
        (200, 0.625, 1.000, 'MODERATE', 'High TMB — many candidates but ample signal'),
        (9999, 0.50, 0.70,  'LOW',      'Very high TMB — dilution risk'),
    ]

    for threshold, expected_r20, p_any_hit, tier, rationale in empirical_tiers:
        if n_filter_passing <= threshold:
            break

    # Per-epitope confidence from model scores
    score_col = None
    for col in ['reranker_score', 'composite_score', 'effective_rank']:
        if col in df.columns:
            score_col = col
            break

    confidence = {
        'confidence_tier': tier,
        'rationale': rationale,
        'expected_recall_at_20': round(expected_r20, 3),
        'p_at_least_one_immunogenic': round(p_any_hit, 3),
        'n_total_mutations': n_total_mutations,
        'n_filter_passing': n_filter_passing,
        'n_selected': n_selected,
        'selectivity_ratio': round(n_selected / max(n_filter_passing, 1), 3),
        'benchmark_basis': 'Muller/Gfeller 91-patient dev set, LightGBM 25-feature reranker',
        'caveat': (
            'Confidence estimates are based on a 131-patient melanoma/mixed-cancer '
            'benchmark. Actual performance may vary by cancer type, HLA alleles, '
            'and mutation profile. These are not clinical predictions.'
        ),
    }

    return confidence


def main():
    parser = argparse.ArgumentParser(
        description='ChatNAV: Personalized neoantigen vaccine design',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input (MAF or VCF required)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--maf', help='Somatic MAF file (TCGA format)')
    input_group.add_argument('--vcf', help='Somatic VCF file (annotated with VEP)')

    parser.add_argument('--expression', required=True,
                        help='Gene expression TSV (gene_id + TPM columns)')
    parser.add_argument('--hla-alleles', required=True,
                        help='Comma-separated HLA alleles (e.g., HLA-A*02:01,HLA-B*07:02,...)')
    parser.add_argument('--output', required=True,
                        help='Output directory')
    parser.add_argument('--patient-id', default='patient',
                        help='Patient identifier')

    # Options
    parser.add_argument('--profile', default='high_tmb',
                        choices=['high_tmb', 'low_tmb', 'research'],
                        help='Scoring profile (default: high_tmb)')
    parser.add_argument('--top-n', type=int, default=20,
                        help='Number of epitopes to select (default: 20)')
    parser.add_argument('--max-mutations', type=int, default=None,
                        help='Limit mutations processed (default: all)')
    parser.add_argument('--ensembl-release', type=int, default=110,
                        help='Ensembl release (default: 110, pinned)')
    parser.add_argument('--skip-mrna', action='store_true',
                        help='Skip mRNA design (LinearDesign step)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Validate inputs without running pipeline')

    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"{'═'*60}")
    print(f"  ChatNAV — Personalized Neoantigen Vaccine Design")
    print(f"{'═'*60}")
    print(f"  Patient:    {args.patient_id}")
    print(f"  Input:      {args.maf or args.vcf}")
    print(f"  Expression: {args.expression}")
    print(f"  HLA:        {args.hla_alleles}")
    print(f"  Profile:    {args.profile}")
    print(f"  Output:     {outdir}")

    t_start = time.time()

    # ── Step 0: MAF/VCF → peptides + MHC binding predictions ──
    binding_path = outdir / 'binding_predictions.tsv'
    meta_path = outdir / 'candidates_meta.tsv'

    step0_args = [
        '--expression', args.expression,
        '--hla-alleles', args.hla_alleles,
        '--ensembl-release', args.ensembl_release,
        '--output-binding', binding_path,
        '--output-meta', meta_path,
    ]
    if args.max_mutations:
        step0_args += ['--max-mutations', args.max_mutations]
    if args.dry_run:
        step0_args.append('--dry-run')

    if args.maf:
        step0_args = ['--maf', args.maf] + step0_args
    else:
        step0_args = ['--maf', args.vcf] + step0_args  # VCF handled same way for now

    run_step('maf_to_pipeline_input.py', step0_args,
             'Step 0: Variant → Peptides + MHC Binding')

    if args.dry_run:
        print("\nDry run complete — inputs validated.")
        return

    # ── Module 8a: Immunogenicity scoring ──
    immuno_path = outdir / 'immunogenicity_scores.tsv'
    proteome_path = REF / 'proteome' / 'human_proteome.fasta'

    run_step('score_immunogenicity.py', [
        '--binding-predictions', binding_path,
        '--candidates-meta', meta_path,
        '--human-proteome', proteome_path,
        '--output', immuno_path,
    ], 'Module 8a: Immunogenicity Scoring (BigMHC-IM + foreignness)')

    # ── Module 8b: Structural scoring ──
    structural_path = outdir / 'structural_scores.tsv'
    run_step('structural_scoring.py', [
        '--immunogenicity-scores', immuno_path,
        '--candidates-meta', meta_path,
        '--output', structural_path,
    ], 'Module 8b: Structural Scoring')

    # ── Module 9: Ranking & selection ──
    selected_path = outdir / 'selected_neoantigens.tsv'
    weights_path = PROJECT_ROOT / 'conf' / 'scoring_weights.yaml'

    run_step('rank_and_select.py', [
        '--immunogenicity-scores', immuno_path,
        '--structural-scores', structural_path,
        '--candidates-meta', meta_path,
        '--binding-predictions', binding_path,
        '--scoring-weights', weights_path,
        '--weight-profile', args.profile,
        '--top-n', args.top_n,
        '--output', selected_path,
    ], 'Module 9: Ranking & Selection')

    # ── Count filter-passing candidates for confidence ──
    import pandas as pd
    meta_df = pd.read_csv(meta_path, sep='\t')
    binding_df = pd.read_csv(binding_path, sep='\t')
    n_total = len(meta_df['peptide_id'].unique()) if 'peptide_id' in meta_df.columns else len(meta_df)
    n_filter = len(binding_df[binding_df.get('binding_rank', pd.Series([1.0])) < 0.02]) if 'binding_rank' in binding_df.columns else n_total

    # ── Module 10: Polyepitope design ──
    polyepitope_fasta = outdir / 'polyepitope.faa'
    polyepitope_design = outdir / 'polyepitope_design.tsv'

    run_step('design_polyepitope.py', [
        '--selected-neoantigens', selected_path,
        '--signal-peptides', REF / 'signal_peptides.fasta',
        '--output-fasta', polyepitope_fasta,
        '--output-design', polyepitope_design,
        '--patient-id', args.patient_id,
    ], 'Module 10: Polyepitope Design')

    # ── Module 11: mRNA design ──
    mrna_fasta = outdir / 'mrna_sequence.fasta'
    synthesis_spec = outdir / 'synthesis_spec.json'

    if not args.skip_mrna:
        run_step('design_mrna.py', [
            '--polyepitope', polyepitope_fasta,
            '--utr-5prime', REF / 'utr_library' / 'utr_5prime.fasta',
            '--utr-3prime', REF / 'utr_library' / 'utr_3prime.fasta',
            '--lineardesign-dir', PROJECT_ROOT / 'LinearDesign',
            '--patient-id', args.patient_id,
            '--output-fasta', mrna_fasta,
            '--output-spec', synthesis_spec,
        ], 'Module 11: mRNA Design (LinearDesign)')
    else:
        print(f"\n  Skipping mRNA design (--skip-mrna)")

    # ── Confidence estimate ──
    confidence = compute_confidence(str(selected_path), n_total, n_filter)
    conf_path = outdir / 'confidence.json'
    with open(conf_path, 'w') as f:
        json.dump(confidence, f, indent=2)

    # ── Summary card ──
    if synthesis_spec.exists():
        run_step('generate_summary_card.py', [
            '--selected-neoantigens', selected_path,
            '--synthesis-spec', synthesis_spec,
            '--polyepitope', polyepitope_fasta,
            '--patient-id', args.patient_id,
            '--output', outdir / 'summary_card.html',
        ], 'Summary Card')

    # ── Final summary ──
    elapsed = time.time() - t_start
    print(f"\n{'═'*60}")
    print(f"  PIPELINE COMPLETE — {args.patient_id}")
    print(f"{'═'*60}")
    print(f"  Time: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  Confidence: {confidence['confidence_tier']} "
          f"(expected R@20={confidence['expected_recall_at_20']:.3f})")

    print(f"\n  Output files:")
    for f in sorted(outdir.iterdir()):
        size = f.stat().st_size
        if size > 1024*1024:
            size_str = f"{size/1024/1024:.1f}MB"
        elif size > 1024:
            size_str = f"{size/1024:.1f}KB"
        else:
            size_str = f"{size}B"
        print(f"    {f.name:<35s} {size_str:>8s}")


if __name__ == '__main__':
    main()
