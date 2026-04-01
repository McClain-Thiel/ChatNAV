#!/usr/bin/env bash
set -euo pipefail

# ══════════════════════════════════════════════════════════════════
# Neoantigen Vaccine Pipeline — Demo Data Preparation
#
# Downloads open-access TCGA-SKCM data for 4 melanoma demo patients.
# All data is open-access (no dbGaP approval required).
#
# Sources:
#   - Masked Somatic Mutation MAFs: GDC Data Portal (open access)
#   - Gene Expression Quantification: GDC Data Portal (STAR counts → TPM)
#   - HLA types: GDC Pan-Immune OptiType calls (open access)
#   - MC3 unified somatic mutations: GDC MC3 publication (open access)
#
# Usage:
#   cd cancer-vaccines
#   bash demo/prepare_demo.sh
# ══════════════════════════════════════════════════════════════════

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$DEMO_DIR")"

# Demo patients: high-TMB TCGA-SKCM melanoma cases
PATIENTS=("TCGA-ER-A19O" "TCGA-EE-A3J5" "TCGA-GN-A262" "TCGA-D9-A4Z6")

GDC_API="https://api.gdc.cancer.gov"

echo "═══════════════════════════════════════════════════"
echo "  Neoantigen Vaccine Pipeline — Demo Data Setup"
echo "═══════════════════════════════════════════════════"

# ── Step 0: Check dependencies ──
for cmd in curl python3 jq; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd is required but not found"
        exit 1
    fi
done

# ── Step 1: Download MC3 Pan-Cancer Somatic Mutations ──
MC3_FILE="$DEMO_DIR/mc3.v0.2.8.PUBLIC.maf.gz"
if [ ! -f "$MC3_FILE" ]; then
    echo ""
    echo "[1/5] Downloading MC3 pan-cancer somatic mutations (~400MB)..."
    echo "      Source: GDC MC3 publication (open access)"
    curl -L -o "$MC3_FILE" \
        "${GDC_API}/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"
    echo "      Done: $MC3_FILE"
else
    echo "[1/5] MC3 file already exists, skipping download"
fi

# ── Step 2: Download HLA type calls (Pan-Immune OptiType) ──
HLA_FILE="$DEMO_DIR/OptiTypeCallsHLA_20171207.tsv"
if [ ! -f "$HLA_FILE" ]; then
    echo ""
    echo "[2/5] Downloading Pan-Immune HLA type calls..."
    echo "      Source: Thorsson et al. 2018 (Immunity) supplementary data"
    # This file is from the GDC Pan-Immune publication
    # Direct GDC API download of the supplementary file
    curl -L -o "$HLA_FILE" \
        "${GDC_API}/data/d0e72ad8-8878-4a43-8946-5be2aedaae4f" 2>/dev/null || \
    echo "      NOTE: HLA file may need manual download from GDC Pan-Immune page"
    echo "      Done: $HLA_FILE"
else
    echo "[2/5] HLA file already exists, skipping download"
fi

# ── Step 3: Download gene expression for each patient ──
echo ""
echo "[3/5] Downloading gene expression quantification per patient..."

for PATIENT in "${PATIENTS[@]}"; do
    PATIENT_DIR="$DEMO_DIR/tcga_skcm/$PATIENT"
    mkdir -p "$PATIENT_DIR"

    EXPR_FILE="$PATIENT_DIR/expression_raw.tsv"
    if [ ! -f "$EXPR_FILE" ]; then
        echo "      Querying GDC for $PATIENT expression data..."

        # Query GDC API for STAR gene counts
        FILE_ID=$(curl -s "${GDC_API}/files" \
            -d '{
                "filters": {
                    "op": "and",
                    "content": [
                        {"op": "in", "content": {"field": "cases.submitter_id", "value": ["'"$PATIENT"'"]}},
                        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
                        {"op": "in", "content": {"field": "data_format", "value": ["TSV"]}},
                        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-SKCM"]}}
                    ]
                },
                "fields": "file_id",
                "size": 1
            }' \
            -H "Content-Type: application/json" | jq -r '.data.hits[0].file_id // empty')

        if [ -n "$FILE_ID" ]; then
            echo "      Downloading $PATIENT expression ($FILE_ID)..."
            curl -s -o "$EXPR_FILE" "${GDC_API}/data/$FILE_ID"
            echo "      Done: $EXPR_FILE"
        else
            echo "      WARNING: No expression file found for $PATIENT"
        fi
    else
        echo "      $PATIENT expression already exists, skipping"
    fi
done

# ── Step 4: Download masked somatic mutation MAFs per patient ──
echo ""
echo "[4/5] Downloading masked somatic mutation MAFs per patient..."

for PATIENT in "${PATIENTS[@]}"; do
    PATIENT_DIR="$DEMO_DIR/tcga_skcm/$PATIENT"
    mkdir -p "$PATIENT_DIR"

    MAF_FILE="$PATIENT_DIR/somatic_masked.maf.gz"
    if [ ! -f "$MAF_FILE" ]; then
        echo "      Querying GDC for $PATIENT masked somatic mutations..."

        FILE_ID=$(curl -s "${GDC_API}/files" \
            -d '{
                "filters": {
                    "op": "and",
                    "content": [
                        {"op": "in", "content": {"field": "cases.submitter_id", "value": ["'"$PATIENT"'"]}},
                        {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
                        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-SKCM"]}}
                    ]
                },
                "fields": "file_id,file_size",
                "size": 1
            }' \
            -H "Content-Type: application/json" | jq -r '.data.hits[0].file_id // empty')

        if [ -n "$FILE_ID" ]; then
            echo "      Downloading $PATIENT MAF ($FILE_ID)..."
            curl -s -o "$MAF_FILE" "${GDC_API}/data/$FILE_ID"
            echo "      Done: $MAF_FILE"
        else
            echo "      WARNING: No MAF file found for $PATIENT, will extract from MC3"
        fi
    else
        echo "      $PATIENT MAF already exists, skipping"
    fi
done

# ── Step 5: Process downloaded data into pipeline format ──
echo ""
echo "[5/5] Converting downloaded data to pipeline input format..."

python3 << 'PYEOF'
import gzip
import json
import os
import sys

DEMO_DIR = os.path.dirname(os.path.abspath("__file__"))
if "demo" not in DEMO_DIR:
    DEMO_DIR = os.path.join(os.getcwd(), "demo")

PATIENTS = ["TCGA-ER-A19O", "TCGA-EE-A3J5", "TCGA-GN-A262", "TCGA-D9-A4Z6"]

# ── Convert MC3 MAF to per-patient VCFs ──
mc3_file = os.path.join(DEMO_DIR, "mc3.v0.2.8.PUBLIC.maf.gz")
if os.path.exists(mc3_file):
    print("  Processing MC3 MAF → per-patient VCF...")
    patient_variants = {p: [] for p in PATIENTS}

    with gzip.open(mc3_file, 'rt') as f:
        header = None
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('Hugo_Symbol'):
                header = line.strip().split('\t')
                continue
            if header is None:
                continue

            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue

            row = dict(zip(header, fields))
            barcode = row.get('Tumor_Sample_Barcode', '')

            # Match patient ID (TCGA barcodes are like TCGA-ER-A19O-01A-11D-...)
            for patient in PATIENTS:
                if barcode.startswith(patient):
                    patient_variants[patient].append(row)
                    break

    for patient, variants in patient_variants.items():
        if not variants:
            print(f"  WARNING: No MC3 variants found for {patient}")
            continue

        patient_dir = os.path.join(DEMO_DIR, "tcga_skcm", patient)
        os.makedirs(patient_dir, exist_ok=True)
        vcf_path = os.path.join(patient_dir, "somatic_annotated.vcf.gz")

        # Write VCF
        with gzip.open(vcf_path, 'wt') as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##source=MC3_v0.2.8\n")
            vcf.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP consequence">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            vcf.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
            vcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
            vcf.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
            vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{patient}_TUMOR\n")

            for v in variants:
                chrom = v.get('Chromosome', '')
                pos = v.get('Start_Position', '')
                ref = v.get('Reference_Allele', '')
                alt = v.get('Tumor_Seq_Allele2', '')
                gene = v.get('Hugo_Symbol', '')
                consequence = v.get('Variant_Classification', '')
                protein = v.get('HGVSp_Short', '').replace('p.', '')

                # Build CSQ-like annotation
                csq = f"{consequence}|{consequence}||{gene}||||||||{protein}"

                # Allele counts
                t_ref = v.get('t_ref_count', '30')
                t_alt = v.get('t_alt_count', '10')
                try:
                    dp = int(t_ref) + int(t_alt)
                    af = int(t_alt) / max(dp, 1)
                except (ValueError, ZeroDivisionError):
                    dp, af = 40, 0.25

                info = f"CSQ={csq}"
                fmt = "GT:AD:DP:AF"
                sample = f"0/1:{t_ref},{t_alt}:{dp}:{af:.3f}"

                vcf.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\t{fmt}\t{sample}\n")

        print(f"  {patient}: {len(variants)} somatic variants → {vcf_path}")

else:
    print("  WARNING: MC3 file not found — VCFs will need to be downloaded separately")

# ── Extract HLA alleles per patient ──
hla_file = os.path.join(DEMO_DIR, "OptiTypeCallsHLA_20171207.tsv")
if os.path.exists(hla_file):
    print("  Processing OptiType HLA calls → per-patient allele files...")
    import csv

    with open(hla_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            barcode = row.get('sample', row.get('Sample', ''))
            for patient in PATIENTS:
                if barcode.startswith(patient):
                    patient_dir = os.path.join(DEMO_DIR, "tcga_skcm", patient)
                    os.makedirs(patient_dir, exist_ok=True)
                    hla_out = os.path.join(patient_dir, "hla_alleles.txt")

                    alleles = []
                    for col in ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']:
                        val = row.get(col, '')
                        if val and val != 'NA':
                            # Normalize to HLA-X*XX:XX format
                            if not val.startswith('HLA-'):
                                gene = col[0]
                                val = f"HLA-{gene}*{val}"
                            alleles.append(val)

                    if alleles:
                        with open(hla_out, 'w') as out:
                            for a in sorted(set(alleles)):
                                out.write(a + '\n')
                        print(f"  {patient}: {len(alleles)} HLA alleles → {hla_out}")
                    break
else:
    print("  WARNING: HLA file not found — creating default HLA alleles")
    # Fall back to common alleles
    for patient in PATIENTS:
        patient_dir = os.path.join(DEMO_DIR, "tcga_skcm", patient)
        os.makedirs(patient_dir, exist_ok=True)
        hla_out = os.path.join(patient_dir, "hla_alleles.txt")
        if not os.path.exists(hla_out):
            with open(hla_out, 'w') as f:
                # Common caucasian HLA alleles as placeholder
                f.write("HLA-A*02:01\n")
                f.write("HLA-A*24:02\n")
                f.write("HLA-B*07:02\n")
                f.write("HLA-B*44:02\n")
                f.write("HLA-C*05:01\n")
                f.write("HLA-C*07:02\n")
            print(f"  {patient}: default HLA alleles (placeholder) → {hla_out}")

# ── Convert expression data to pipeline format ──
print("  Processing expression data → TPM format...")
for patient in PATIENTS:
    patient_dir = os.path.join(DEMO_DIR, "tcga_skcm", patient)
    raw_expr = os.path.join(patient_dir, "expression_raw.tsv")
    tpm_out = os.path.join(patient_dir, "expression_tpm.tsv")

    if os.path.exists(raw_expr) and not os.path.exists(tpm_out):
        try:
            import pandas as pd
            df = pd.read_csv(raw_expr, sep='\t', comment='#')
            # GDC STAR counts file has columns: gene_id, gene_name, gene_type,
            # unstranded, stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded
            if 'tpm_unstranded' in df.columns:
                result = pd.DataFrame({
                    'gene_id': df['gene_id'],
                    'gene_name': df['gene_name'],
                    'tpm': df['tpm_unstranded'],
                    'expression_imputed': False,
                })
            elif 'TPM' in df.columns:
                result = df[['gene_id', 'gene_name', 'TPM']].copy()
                result.columns = ['gene_id', 'gene_name', 'tpm']
                result['expression_imputed'] = False
            else:
                # Try to compute TPM from counts
                print(f"  WARNING: No TPM column in {raw_expr}, using raw counts")
                continue

            # Filter out non-gene rows (metadata rows at top of GDC files)
            result = result[result['gene_id'].str.startswith('ENSG', na=False)]
            result.to_csv(tpm_out, sep='\t', index=False)
            print(f"  {patient}: {len(result)} genes with TPM → {tpm_out}")
        except Exception as e:
            print(f"  WARNING: Failed to process {patient} expression: {e}")

print("")
print("═══════════════════════════════════════════════════")
print("  Demo data preparation complete!")
print("")
print("  Next steps:")
print("    1. Start Docker Desktop")
print("    2. docker build -t neoantigen/scoring:0.1.0 docker/scoring/")
print("    3. docker build -t neoantigen/design:0.1.0 docker/design/")
print("    4. nextflow run main.nf -profile local,test")
print("═══════════════════════════════════════════════════")
PYEOF

echo ""
echo "Done!"
