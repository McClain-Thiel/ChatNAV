#!/usr/bin/env python3
"""
Convert a VEP-annotated somatic VCF to MAF format for the ChatNAV pipeline.

Supports:
  - Single-sample VCF (tumor-only)
  - Paired VCF (tumor + normal, e.g., from Mutect2/Strelka2)
  - VEP CSQ annotations (required for protein change)
  - Gzipped or plain VCF

Extracts the columns needed by maf_to_pipeline_input.py:
  Hugo_Symbol, Chromosome, Start_Position, End_Position,
  Variant_Classification, HGVSp_Short, Tumor_Seq_Allele2,
  Reference_Allele, t_alt_count, t_ref_count

Usage:
    python bin/vcf_to_maf.py \
        --vcf patient_somatic.vcf.gz \
        --output patient_somatic.maf \
        --tumor-sample TUMOR \
        --normal-sample NORMAL
"""

import argparse
import gzip
import re
import sys
from pathlib import Path


# VEP Variant_Classification mapping
VEP_TO_MAF_CLASS = {
    'missense_variant': 'Missense_Mutation',
    'frameshift_variant': 'Frame_Shift_Del',
    'inframe_deletion': 'In_Frame_Del',
    'inframe_insertion': 'In_Frame_Ins',
    'stop_gained': 'Nonsense_Mutation',
    'stop_lost': 'Nonstop_Mutation',
    'start_lost': 'Translation_Start_Site',
    'splice_acceptor_variant': 'Splice_Site',
    'splice_donor_variant': 'Splice_Site',
    'synonymous_variant': 'Silent',
    'intron_variant': 'Intron',
    '3_prime_utr_variant': "3'UTR",
    '5_prime_utr_variant': "5'UTR",
    'upstream_gene_variant': "5'Flank",
    'downstream_gene_variant': "3'Flank",
    'intergenic_variant': 'IGR',
}


def parse_csq_header(header_line: str) -> list[str]:
    """Extract CSQ field names from VCF header line.

    Looks for: ##INFO=<ID=CSQ,...,Description="...Format: Allele|Consequence|...">
    """
    match = re.search(r'Format:\s*([^"]+)', header_line)
    if not match:
        raise ValueError("Cannot parse CSQ Format from VCF header. "
                         "Is this VCF annotated with VEP?")
    return match.group(1).strip().split('|')


def parse_vcf(vcf_path: str, tumor_sample: str = None, normal_sample: str = None):
    """Parse VCF and yield MAF-format rows."""
    opener = gzip.open if str(vcf_path).endswith('.gz') else open

    csq_fields = None
    sample_names = []
    tumor_idx = None
    normal_idx = None

    with opener(vcf_path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')

            # Header lines
            if line.startswith('##'):
                if 'ID=CSQ' in line:
                    csq_fields = parse_csq_header(line)
                continue

            # Column header
            if line.startswith('#CHROM'):
                cols = line.split('\t')
                sample_names = cols[9:] if len(cols) > 9 else []

                # Find tumor/normal sample indices
                if tumor_sample and tumor_sample in sample_names:
                    tumor_idx = sample_names.index(tumor_sample) + 9
                elif sample_names:
                    tumor_idx = 9  # default to first sample
                if normal_sample and normal_sample in sample_names:
                    normal_idx = sample_names.index(normal_sample) + 9
                continue

            # Data lines
            fields = line.split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0].replace('chr', '')
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(',')
            info = fields[7]
            filter_val = fields[6]

            # Skip filtered variants (keep PASS and .)
            if filter_val not in ('PASS', '.'):
                continue

            # Parse allele depths from tumor sample
            t_ref_count = 0
            t_alt_count = 0
            if tumor_idx and tumor_idx < len(fields):
                fmt = fields[8].split(':')
                sample_data = fields[tumor_idx].split(':')
                if 'AD' in fmt:
                    ad_idx = fmt.index('AD')
                    if ad_idx < len(sample_data):
                        ad = sample_data[ad_idx].split(',')
                        t_ref_count = int(ad[0]) if ad[0] != '.' else 0
                        t_alt_count = int(ad[1]) if len(ad) > 1 and ad[1] != '.' else 0

            # Parse CSQ annotations
            csq_entries = []
            for part in info.split(';'):
                if part.startswith('CSQ='):
                    for entry in part[4:].split(','):
                        csq_entries.append(entry.split('|'))

            if not csq_entries or not csq_fields:
                continue

            # Pick the most severe consequence per alt allele
            for alt in alts:
                best_entry = None
                best_rank = 999

                for entry in csq_entries:
                    if len(entry) != len(csq_fields):
                        continue
                    csq = dict(zip(csq_fields, entry))

                    # Match allele
                    if csq.get('Allele', '') != alt and len(alts) > 1:
                        continue

                    # Rank by consequence severity
                    consequence = csq.get('Consequence', '')
                    for i, (vep_term, _) in enumerate(VEP_TO_MAF_CLASS.items()):
                        if vep_term in consequence:
                            if i < best_rank:
                                best_rank = i
                                best_entry = csq
                            break

                if best_entry is None:
                    continue

                consequence = best_entry.get('Consequence', '')
                maf_class = 'Unknown'
                for vep_term, maf_term in VEP_TO_MAF_CLASS.items():
                    if vep_term in consequence:
                        maf_class = maf_term
                        break

                # Determine frameshift direction
                if maf_class == 'Frame_Shift_Del' and len(alt) > len(ref):
                    maf_class = 'Frame_Shift_Ins'

                # Extract HGVSp
                hgvsp = best_entry.get('HGVSp', '')
                if ':' in hgvsp:
                    hgvsp = hgvsp.split(':')[1]

                # Convert 3-letter to 1-letter amino acid codes
                hgvsp_short = convert_hgvsp_to_short(hgvsp)

                gene = best_entry.get('SYMBOL', best_entry.get('Gene', ''))

                # Compute end position
                end_pos = pos + len(ref) - 1

                yield {
                    'Hugo_Symbol': gene,
                    'Chromosome': chrom,
                    'Start_Position': pos,
                    'End_Position': end_pos,
                    'Reference_Allele': ref,
                    'Tumor_Seq_Allele2': alt,
                    'Variant_Classification': maf_class,
                    'HGVSp_Short': hgvsp_short,
                    'Variant_Type': classify_variant_type(ref, alt),
                    't_alt_count': t_alt_count,
                    't_ref_count': t_ref_count,
                }


AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*',
}


def convert_hgvsp_to_short(hgvsp: str) -> str:
    """Convert VEP HGVSp (3-letter) to short (1-letter) format.

    e.g., p.Ser592Phe → p.S592F
          p.Tyr490LeufsTer4 → p.Y490Lfs*4
    """
    if not hgvsp or not hgvsp.startswith('p.'):
        return hgvsp

    result = 'p.'
    remaining = hgvsp[2:]

    # Replace 3-letter codes with 1-letter
    i = 0
    while i < len(remaining):
        matched = False
        for aa3, aa1 in AA_3TO1.items():
            if remaining[i:i+3] == aa3:
                result += aa1
                i += 3
                matched = True
                break
        if not matched:
            result += remaining[i]
            i += 1

    return result


def classify_variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    elif len(ref) > len(alt):
        return 'DEL'
    elif len(ref) < len(alt):
        return 'INS'
    else:
        return 'DNP' if len(ref) == 2 else 'ONP'


def main():
    parser = argparse.ArgumentParser(description='Convert VEP-annotated VCF to MAF')
    parser.add_argument('--vcf', required=True, help='Input VCF (VEP-annotated)')
    parser.add_argument('--output', required=True, help='Output MAF file')
    parser.add_argument('--tumor-sample', default=None,
                        help='Tumor sample name in VCF (default: first sample)')
    parser.add_argument('--normal-sample', default=None,
                        help='Normal sample name in VCF (for paired calling)')
    args = parser.parse_args()

    if not Path(args.vcf).exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)

    rows = list(parse_vcf(args.vcf, args.tumor_sample, args.normal_sample))

    if not rows:
        print("ERROR: No variants parsed from VCF. Is it VEP-annotated (CSQ field)?",
              file=sys.stderr)
        sys.exit(1)

    # Write MAF
    maf_cols = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification',
        'Variant_Type', 'HGVSp_Short', 't_alt_count', 't_ref_count',
    ]

    with open(args.output, 'w') as f:
        f.write('\t'.join(maf_cols) + '\n')
        for row in rows:
            f.write('\t'.join(str(row.get(c, '')) for c in maf_cols) + '\n')

    # Summary
    from collections import Counter
    class_counts = Counter(r['Variant_Classification'] for r in rows)
    print(f"Converted {len(rows)} variants from VCF → MAF", file=sys.stderr)
    for cls, count in class_counts.most_common():
        print(f"  {cls}: {count}", file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
