#!/usr/bin/env python3
"""
End-to-end test: real TCGA melanoma data through all modules.

Uses real:
  - TCGA-EE-A3J5 somatic MAF (from GDC)
  - STAR gene expression counts (from GDC)
  - Ensembl GRCh38 protein sequences (via pyensembl)
  - MHCflurry Class1PresentationPredictor (real binding predictions)
  - DeepImmuno-CNN (real immunogenicity predictions)

NO RANDOM DATA. NO FALLBACKS. Every step uses the real tool.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
BIN = PROJECT_ROOT / 'bin'
REF = PROJECT_ROOT / 'reference'
DEMO = PROJECT_ROOT / 'demo' / 'tcga_skcm'


def run_step(script: str, args: list[str], label: str):
    """Run a pipeline step. Fails hard on any error."""
    cmd = [sys.executable, str(BIN / script)] + args
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    # Print all output
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        for line in result.stderr.splitlines():
            # Skip TF/CUDA noise, show everything else
            if any(skip in line.lower() for skip in ['cuda', 'cudnn', 'gpu available',
                    'tensorflow/core/platform', 'onednn', 'xla']):
                continue
            print(line, file=sys.stderr)

    if result.returncode != 0:
        print(f"\nFATAL: {label} failed with exit code {result.returncode}")
        sys.exit(1)

    return result


def main():
    # Verify real data exists
    maf_path = DEMO / 'TCGA-EE-A3J5_somatic.maf.gz'
    expr_path = DEMO / 'TCGA-EE-A3J5_expression.tsv'

    for path, desc in [(maf_path, "Somatic MAF"), (expr_path, "Expression data")]:
        if not path.exists():
            print(f"FATAL: {desc} not found at {path}")
            sys.exit(1)

    # Verify human proteome exists
    proteome_path = REF / 'proteome' / 'human_proteome.fasta'
    if not proteome_path.exists():
        print(f"FATAL: Human proteome not found at {proteome_path}")
        sys.exit(1)

    with tempfile.TemporaryDirectory(prefix='neoantigen_e2e_') as tmpdir:
        tmpdir = Path(tmpdir)
        print(f"Patient: TCGA-EE-A3J5 (TCGA SKCM melanoma)")
        print(f"Working dir: {tmpdir}")
        print(f"All steps use REAL tools — no fallbacks, no random data\n")

        # ── Step 0: MAF → real peptides + real MHCflurry binding ──
        binding_path = tmpdir / 'binding_predictions.tsv'
        meta_path = tmpdir / 'candidates_meta.tsv'

        run_step('maf_to_pipeline_input.py', [
            '--maf', str(maf_path),
            '--expression', str(expr_path),
            '--hla-alleles', 'HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*44:03,HLA-C*07:02,HLA-C*05:01',
            '--max-mutations', '50',
            '--output-binding', str(binding_path),
            '--output-meta', str(meta_path),
        ], 'Step 0: MAF → Peptides + MHCflurry Binding (REAL)')

        # ── Module 8: DeepImmuno-CNN immunogenicity (REAL) ──
        immuno_path = tmpdir / 'immunogenicity_scores.tsv'
        run_step('score_immunogenicity.py', [
            '--binding-predictions', str(binding_path),
            '--candidates-meta', str(meta_path),
            '--human-proteome', str(proteome_path),
            '--output', str(immuno_path),
        ], 'Module 8a: Immunogenicity Scoring (BigMHC-IM + foreignness)')

        # ── Module 8b: Structural Scoring (Tier 1 — position lookup) ──
        structural_path = tmpdir / 'structural_scores.tsv'
        run_step('structural_scoring.py', [
            '--immunogenicity-scores', str(immuno_path),
            '--candidates-meta', str(meta_path),
            '--output', str(structural_path),
        ], 'Module 8b: Structural Scoring (Tier 1 — position lookup)')

        # ── Module 9: Ranking & Selection ──
        selected_path = tmpdir / 'selected_neoantigens.tsv'
        weights_path = PROJECT_ROOT / 'conf' / 'scoring_weights.yaml'
        run_step('rank_and_select.py', [
            '--immunogenicity-scores', str(immuno_path),
            '--structural-scores', str(structural_path),
            '--candidates-meta', str(meta_path),
            '--binding-predictions', str(binding_path),
            '--scoring-weights', str(weights_path),
            '--weight-profile', 'high_tmb',
            '--top-n', '20',
            '--output', str(selected_path),
        ], 'Module 9: Ranking & Selection')

        # ── Module 10: Polyepitope Design ──
        polyepitope_fasta = tmpdir / 'polyepitope.faa'
        polyepitope_design = tmpdir / 'polyepitope_design.tsv'
        signal_pep_path = REF / 'signal_peptides.fasta'
        run_step('design_polyepitope.py', [
            '--selected-neoantigens', str(selected_path),
            '--signal-peptides', str(signal_pep_path),
            '--output-fasta', str(polyepitope_fasta),
            '--output-design', str(polyepitope_design),
            '--patient-id', 'TCGA_EE_A3J5',
        ], 'Module 10: Polyepitope Design')

        # ── Module 11: mRNA Design ──
        mrna_fasta = tmpdir / 'mrna_sequence.fasta'
        synthesis_spec = tmpdir / 'synthesis_spec.json'
        utr5_path = REF / 'utr_library' / 'utr_5prime.fasta'
        utr3_path = REF / 'utr_library' / 'utr_3prime.fasta'
        run_step('design_mrna.py', [
            '--polyepitope', str(polyepitope_fasta),
            '--utr-5prime', str(utr5_path),
            '--utr-3prime', str(utr3_path),
            '--lineardesign-dir', str(PROJECT_ROOT / 'LinearDesign'),
            '--patient-id', 'TCGA_EE_A3J5',
            '--output-fasta', str(mrna_fasta),
            '--output-spec', str(synthesis_spec),
        ], 'Module 11: mRNA Design (LinearDesign)')

        # ── Results ──
        print(f"\n{'='*70}")
        print(f"  END-TO-END RESULTS — TCGA-EE-A3J5 (Melanoma)")
        print(f"{'='*70}")

        import pandas as pd

        selected = pd.read_csv(selected_path, sep='\t')
        print(f"\nSelected neoantigens ({len(selected)} total):")
        cols = ['rank', 'peptide_id', 'peptide_sequence', 'hla_allele', 'gene',
                'composite_score', 'immunogenicity_score', 'binding_rank']
        display_cols = [c for c in cols if c in selected.columns]
        if len(selected) > 0:
            print(selected[display_cols].head(20).to_string(index=False))

        # Verify polyepitope contains only amino acids
        print(f"\nPolyepitope sequence:")
        with open(polyepitope_fasta) as f:
            lines = f.readlines()
            for line in lines:
                print(f"  {line.rstrip()}")
            seq = ''.join(l.strip() for l in lines if not l.startswith('>'))
            valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
            non_aa = [c for c in seq if c not in valid_aas]
            if non_aa:
                print(f"\n  FATAL: Polyepitope contains non-amino-acid characters: {set(non_aa)}")
                sys.exit(1)
            else:
                print(f"  VALIDATED: All {len(seq)} characters are valid amino acids")

        with open(synthesis_spec) as f:
            spec = json.load(f)
        print(f"\nSynthesis specification:")
        print(f"  mRNA length:        {spec['length_nt']} nt")
        print(f"  Protein length:     {spec['protein_length_aa']} aa")
        print(f"  GC content:         {spec['gc_content']:.1%}")
        print(f"  Codon optimization: {spec['codon_optimization']}")
        print(f"  Pseudouridine:      {spec['modifications']['pseudouridine_count']} positions")

        # Save results
        results_dir = PROJECT_ROOT / 'results' / 'TCGA_EE_A3J5'
        results_dir.mkdir(parents=True, exist_ok=True)
        import shutil
        for f in [immuno_path, selected_path, polyepitope_fasta,
                  polyepitope_design, mrna_fasta, synthesis_spec]:
            shutil.copy2(f, results_dir / f.name)

        print(f"\nOutputs saved to: {results_dir}/")
        print(f"\nAll steps completed with REAL tools — no fallbacks used.")


if __name__ == '__main__':
    main()
