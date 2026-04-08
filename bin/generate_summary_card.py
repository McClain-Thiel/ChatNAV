#!/usr/bin/env python3
"""
Vaccine Summary Card Generator

Produces a one-page self-contained HTML report for a patient's vaccine design.
Reads pipeline outputs and renders:
  - Patient info + HLA alleles
  - Top 20 epitopes table with scores and flags
  - Polyepitope sequence with linkers highlighted
  - mRNA stats (length, GC%, MFE, CAI)
  - Filter funnel counts
  - CDS quality scan results

Usage:
    python bin/generate_summary_card.py \
        --selected-neoantigens results/patient/selected_neoantigens.tsv \
        --synthesis-spec results/patient/synthesis_spec.json \
        --polyepitope results/patient/polyepitope.faa \
        --patient-id TCGA-EE-A3J5 \
        --output results/patient/summary_card.html
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def render_epitope_table(df: pd.DataFrame) -> str:
    """Render the selected neoantigens as an HTML table."""
    cols = ['rank', 'peptide_sequence', 'hla_allele', 'gene', 'mutation',
            'binding_rank', 'immunogenicity_score', 'validation_flags']
    display_cols = [c for c in cols if c in df.columns]

    rows_html = []
    for _, row in df.iterrows():
        mhc_class = str(row.get('mhc_class', 'I'))
        row_class = 'class-ii' if mhc_class == 'II' else ''
        cells = []
        for col in display_cols:
            val = row.get(col, '')
            if isinstance(val, float):
                val = f'{val:.4f}' if col in ('binding_rank', 'immunogenicity_score') else f'{val:.2f}'
            cells.append(f'<td>{val}</td>')
        rows_html.append(f'<tr class="{row_class}">{"".join(cells)}</tr>')

    header = ''.join(f'<th>{c.replace("_", " ").title()}</th>' for c in display_cols)
    return f'''<table>
<thead><tr>{header}</tr></thead>
<tbody>{"".join(rows_html)}</tbody>
</table>'''


def render_polyepitope(fasta_path: str) -> str:
    """Render polyepitope sequence with linkers highlighted."""
    seq = ''
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()

    # Highlight linkers
    highlighted = seq.replace('AAY', '<span class="linker-i">AAY</span>')
    highlighted = highlighted.replace('GPGPG', '<span class="linker-ii">GPGPG</span>')
    highlighted = highlighted.replace('AAAAA', '<span class="linker-term">AAAAA</span>')

    # Wrap at 60 chars (preserving HTML tags)
    return f'<div class="sequence">{highlighted}</div>'


def render_mrna_stats(spec: dict) -> str:
    """Render mRNA statistics."""
    qc = spec.get('quality_scans', {})
    return f'''<table class="stats">
<tr><td>Total length</td><td>{spec["length_nt"]} nt</td></tr>
<tr><td>Protein length</td><td>{spec["protein_length_aa"]} aa</td></tr>
<tr><td>GC content</td><td>{spec["gc_content"]:.1%}</td></tr>
<tr><td>Codon optimization</td><td>{spec["codon_optimization"]}</td></tr>
<tr><td>MFE</td><td>{spec.get("lineardesign_mfe_kcal", "N/A")} kcal/mol</td></tr>
<tr><td>CAI</td><td>{spec.get("lineardesign_cai", "N/A")}</td></tr>
<tr><td>Pseudouridine positions</td><td>{spec["modifications"]["pseudouridine_count"]}</td></tr>
<tr><td>Cap analog</td><td>{spec["modifications"]["cap_analog"]}</td></tr>
<tr><td>Poly-A tail</td><td>{spec["modifications"]["poly_a_length"]} nt</td></tr>
<tr><td>Internal poly-A signals</td><td>{qc.get("internal_poly_a_signals", "not scanned")}</td></tr>
<tr><td>Cryptic splice sites</td><td>{qc.get("cryptic_splice_sites", "not scanned")}</td></tr>
</table>'''


HTML_TEMPLATE = '''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Vaccine Summary — {patient_id}</title>
<style>
  @media print {{ @page {{ margin: 1cm; }} }}
  body {{ font-family: -apple-system, "Segoe UI", Roboto, sans-serif; max-width: 1100px; margin: 0 auto; padding: 20px; color: #1a1a1a; font-size: 13px; }}
  h1 {{ font-size: 20px; border-bottom: 2px solid #2563eb; padding-bottom: 8px; margin-bottom: 4px; }}
  h2 {{ font-size: 15px; margin: 16px 0 8px; color: #2563eb; }}
  .subtitle {{ color: #666; font-size: 13px; margin-bottom: 16px; }}
  .hla {{ font-family: monospace; background: #f1f5f9; padding: 2px 6px; border-radius: 3px; margin: 2px; display: inline-block; }}
  .hla-ii {{ background: #fef3c7; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 12px; margin: 8px 0; }}
  th {{ background: #f8fafc; text-align: left; padding: 6px 8px; border: 1px solid #e2e8f0; font-weight: 600; }}
  td {{ padding: 5px 8px; border: 1px solid #e2e8f0; }}
  tr:nth-child(even) {{ background: #f8fafc; }}
  tr.class-ii {{ background: #fefce8; }}
  tr.class-ii td {{ border-left-color: #f59e0b; }}
  .stats td:first-child {{ font-weight: 600; width: 200px; }}
  .stats {{ width: auto; }}
  .sequence {{ font-family: "Courier New", monospace; font-size: 12px; word-break: break-all; line-height: 1.8; background: #f8fafc; padding: 12px; border-radius: 6px; border: 1px solid #e2e8f0; }}
  .linker-i {{ background: #dbeafe; padding: 1px 3px; border-radius: 2px; font-weight: bold; }}
  .linker-ii {{ background: #fef3c7; padding: 1px 3px; border-radius: 2px; font-weight: bold; }}
  .linker-term {{ background: #dcfce7; padding: 1px 3px; border-radius: 2px; font-weight: bold; }}
  .legend {{ display: flex; gap: 16px; font-size: 11px; margin: 6px 0; }}
  .legend span {{ padding: 2px 8px; border-radius: 3px; }}
  .cols {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
  .footer {{ margin-top: 20px; padding-top: 10px; border-top: 1px solid #e2e8f0; font-size: 11px; color: #94a3b8; }}
</style>
</head>
<body>
<h1>Vaccine Summary Card</h1>
<div class="subtitle">Patient: <strong>{patient_id}</strong> | Generated: {date}</div>

<h2>HLA Alleles</h2>
<div>{hla_html}</div>

<h2>Selected Neoantigens ({n_epitopes} total: {n_class_i} Class I + {n_class_ii} Class II)</h2>
{epitope_table}

<h2>Polyepitope Construct ({polyepitope_len} aa)</h2>
<div class="legend">
  <span class="linker-i">AAY</span> MHC-I linker
  <span class="linker-ii">GPGPG</span> MHC-II linker
  <span class="linker-term">AAAAA</span> Terminal
</div>
{polyepitope_html}

<div class="cols">
<div>
<h2>mRNA Specification</h2>
{mrna_stats}
</div>
<div>
<h2>Construct Architecture</h2>
<table class="stats">
<tr><td>5' Cap</td><td>CleanCap-AG</td></tr>
<tr><td>5' UTR</td><td>{utr5_len} nt (HBB-derived)</td></tr>
<tr><td>Signal peptide</td><td>21 aa</td></tr>
<tr><td>Epitopes</td><td>{n_epitopes} ({n_class_i} Class I + {n_class_ii} Class II)</td></tr>
<tr><td>MITD domain</td><td>58 aa</td></tr>
<tr><td>3' UTR</td><td>{utr3_len} nt (AES/mtRNR1)</td></tr>
<tr><td>Poly-A tail</td><td>{poly_a_len} nt</td></tr>
</table>
</div>
</div>

<div class="footer">
  Generated by ChatNAV | N1-methylpseudouridine at all U positions | For research use only
</div>
</body>
</html>'''


def main():
    parser = argparse.ArgumentParser(description='Generate vaccine summary card')
    parser.add_argument('--selected-neoantigens', required=True)
    parser.add_argument('--synthesis-spec', required=True)
    parser.add_argument('--polyepitope', required=True)
    parser.add_argument('--patient-id', default='patient')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.selected_neoantigens, sep='\t')
    with open(args.synthesis_spec) as f:
        spec = json.load(f)

    # HLA alleles
    hla_alleles = sorted(df['hla_allele'].dropna().unique()) if 'hla_allele' in df.columns else []
    hla_html = ' '.join(
        f'<span class="hla{"" if not any(x in a for x in ("DRB","DQB","DPB")) else " hla-ii"}">{a}</span>'
        for a in hla_alleles
    )

    # Class counts
    if 'mhc_class' in df.columns:
        n_class_i = len(df[df['mhc_class'] != 'II'])
        n_class_ii = len(df[df['mhc_class'] == 'II'])
    else:
        n_class_i = len(df)
        n_class_ii = 0

    from datetime import date
    html = HTML_TEMPLATE.format(
        patient_id=args.patient_id,
        date=date.today().isoformat(),
        hla_html=hla_html,
        n_epitopes=len(df),
        n_class_i=n_class_i,
        n_class_ii=n_class_ii,
        epitope_table=render_epitope_table(df),
        polyepitope_html=render_polyepitope(args.polyepitope),
        polyepitope_len=spec.get('protein_length_aa', '?'),
        mrna_stats=render_mrna_stats(spec),
        utr5_len=spec.get('utr_5prime_length', '?'),
        utr3_len=spec.get('utr_3prime_length', '?'),
        poly_a_len=spec.get('modifications', {}).get('poly_a_length', '?'),
    )

    with open(args.output, 'w') as f:
        f.write(html)

    print(f"Vaccine summary card → {args.output}")


if __name__ == '__main__':
    main()
