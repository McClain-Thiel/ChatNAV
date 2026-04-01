#!/usr/bin/env python3
"""
Module 5 helper: Prepare PyClone-VI input from annotated VCF + CNV segments.

Extracts per-variant: mutation_id, ref_counts, alt_counts, major_cn, minor_cn, normal_cn
PyClone-VI expects a TSV with these columns.

Inputs:
    somatic_annotated.vcf.gz    — annotated VCF from Module 2
    cnv_segments.tsv            — CNV calls (chr, start, end, log2_ratio, total_cn, minor_cn)

Output:
    pyclone_input.tsv           — PyClone-VI formatted input
"""

import argparse
import gzip
import sys

import pandas as pd


def parse_vcf_variants(vcf_path: str) -> list[dict]:
    """Extract variant info from VCF (supports .gz)."""
    variants = []
    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue

            chrom, pos, vid, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            filt = fields[6]

            # Only PASS variants
            if filt not in ('PASS', '.'):
                continue

            # Parse FORMAT and sample columns
            fmt_keys = fields[8].split(':')
            # Tumor is typically the last sample column
            tumor_values = fields[-1].split(':')
            fmt = dict(zip(fmt_keys, tumor_values))

            # Extract allele depths
            ref_count, alt_count = 0, 0
            if 'AD' in fmt:
                ad = fmt['AD'].split(',')
                ref_count = int(ad[0])
                alt_count = int(ad[1]) if len(ad) > 1 else 0
            elif 'DP' in fmt and 'AF' in fmt:
                dp = int(fmt['DP'])
                af = float(fmt['AF'])
                alt_count = int(dp * af)
                ref_count = dp - alt_count

            # Parse VEP annotation for gene/consequence if present
            info = fields[7]
            gene = ''
            consequence = ''
            protein_change = ''
            for field in info.split(';'):
                if field.startswith('CSQ=') or field.startswith('ANN='):
                    ann = field.split('=')[1].split(',')[0]
                    ann_parts = ann.split('|')
                    if len(ann_parts) > 3:
                        gene = ann_parts[3] if len(ann_parts) > 3 else ''
                        consequence = ann_parts[1] if len(ann_parts) > 1 else ''
                        protein_change = ann_parts[10] if len(ann_parts) > 10 else ''

            mutation_id = f"{chrom}_{pos}_{ref}_{alt}"

            variants.append({
                'mutation_id': mutation_id,
                'chr': chrom,
                'position': pos,
                'ref': ref,
                'alt': alt,
                'ref_counts': ref_count,
                'alt_counts': alt_count,
                'vaf': alt_count / max(ref_count + alt_count, 1),
                'gene': gene,
                'consequence': consequence,
                'protein_change': protein_change,
            })

    return variants


def assign_copy_number(variants: list[dict], cnv_path: str) -> list[dict]:
    """Assign copy number to each variant from CNV segments."""
    if not cnv_path or cnv_path == 'NO_CNV':
        # Default: diploid
        for v in variants:
            v['major_cn'] = 1
            v['minor_cn'] = 1
            v['normal_cn'] = 2
        return variants

    cnv = pd.read_csv(cnv_path, sep='\t')

    # Standardize column names
    cnv.columns = [c.lower().replace(' ', '_') for c in cnv.columns]

    for v in variants:
        # Find overlapping CNV segment
        chrom = v['chr']
        pos = v['position']

        overlap = cnv[
            (cnv['chr'] == chrom) | (cnv['chr'] == chrom.replace('chr', ''))
        ]
        overlap = overlap[(overlap['start'] <= pos) & (overlap['end'] >= pos)]

        if not overlap.empty:
            row = overlap.iloc[0]
            total_cn = int(row.get('total_cn', row.get('ploidy', 2)))
            minor_cn = int(row.get('minor_cn', max(0, total_cn // 2 - 1)))
            major_cn = total_cn - minor_cn
        else:
            major_cn, minor_cn = 1, 1

        v['major_cn'] = major_cn
        v['minor_cn'] = minor_cn
        v['normal_cn'] = 2

    return variants


def main():
    parser = argparse.ArgumentParser(description='Prepare PyClone-VI input')
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--cnv-segments', default=None)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    variants = parse_vcf_variants(args.vcf)
    print(f"Parsed {len(variants)} PASS variants from VCF", file=sys.stderr)

    variants = assign_copy_number(variants, args.cnv_segments)

    # PyClone-VI input format
    rows = []
    for v in variants:
        rows.append({
            'mutation_id': v['mutation_id'],
            'sample_id': 'tumor',
            'ref_counts': v['ref_counts'],
            'alt_counts': v['alt_counts'],
            'major_cn': v['major_cn'],
            'minor_cn': v['minor_cn'],
            'normal_cn': v['normal_cn'],
            # Extra columns for downstream use
            'chr': v['chr'],
            'position': v['position'],
            'ref': v['ref'],
            'alt': v['alt'],
            'vaf': round(v['vaf'], 4),
            'gene': v['gene'],
            'consequence': v['consequence'],
            'protein_change': v['protein_change'],
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(df)} variants to PyClone-VI input: {args.output}")


if __name__ == '__main__':
    main()
