#!/bin/bash
# Restore all pipeline data from S3 bucket
# Run this on a new instance instead of re-downloading everything
set -e

BUCKET="s3://chatnav-pipeline-data"
echo "Restoring pipeline data from $BUCKET"

# MHCflurry models
echo "=== MHCflurry models ==="
MHCFLURRY_DIR="$HOME/.local/share/mhcflurry"
mkdir -p "$MHCFLURRY_DIR"
aws s3 sync "$BUCKET/mhcflurry-models/" "$MHCFLURRY_DIR/" --quiet
echo "Done ($(du -sh $MHCFLURRY_DIR | cut -f1))"

# Ensembl cache
echo "=== Ensembl cache ==="
ENSEMBL_DIR="$HOME/.cache/pyensembl/GRCh38/ensembl110"
mkdir -p "$ENSEMBL_DIR"
aws s3 sync "$BUCKET/ensembl-cache/GRCh38/ensembl110/" "$ENSEMBL_DIR/" --quiet
echo "Done ($(du -sh $ENSEMBL_DIR | cut -f1))"

# Human proteome
echo "=== Human proteome ==="
mkdir -p reference/proteome
aws s3 cp "$BUCKET/reference/human_proteome.fasta" reference/proteome/human_proteome.fasta --quiet
aws s3 cp "$BUCKET/reference/human_proteome.fasta.gz" reference/proteome/human_proteome.fasta.gz --quiet
echo "Done"

# PANDORA database (if available)
echo "=== PANDORA database ==="
PANDORA_DATA=$(python3 -c "import PANDORA; print(PANDORA.PANDORA_data)" 2>/dev/null || echo "")
if [ -n "$PANDORA_DATA" ]; then
    mkdir -p "$PANDORA_DATA"
    aws s3 sync "$BUCKET/pandora-db/" "$PANDORA_DATA/" --quiet
    echo "Done ($(du -sh $PANDORA_DATA | cut -f1))"
else
    echo "PANDORA not installed — skipping"
fi

# Demo data
echo "=== Demo data ==="
mkdir -p demo/tcga_skcm
aws s3 sync "$BUCKET/demo/tcga_skcm/" demo/tcga_skcm/ --quiet
echo "Done"

# Benchmark data
echo "=== Benchmark data ==="
mkdir -p benchmark/gartner
aws s3 cp "$BUCKET/benchmark/gartner/supplementary_table.xlsx" benchmark/gartner/supplementary_table.xlsx --quiet 2>/dev/null || echo "Not available"
echo "Done"

# PANDORA Docker image
echo "=== PANDORA Docker image ==="
if ! docker images pandora | grep -q latest; then
    aws s3 cp "$BUCKET/docker/pandora.tar.gz" /tmp/pandora.tar.gz --quiet
    gunzip -f /tmp/pandora.tar.gz
    docker load < /tmp/pandora.tar
    rm /tmp/pandora.tar
    echo "Loaded"
else
    echo "Already loaded"
fi

echo ""
echo "=== All data restored from S3 ==="
echo "Total local cache size:"
du -sh "$HOME/.local/share/mhcflurry" "$HOME/.cache/pyensembl" reference/proteome 2>/dev/null
