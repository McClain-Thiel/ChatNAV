# Setup Log — Neoantigen Vaccine Pipeline

Everything needed to reproduce this environment from scratch.

## System Requirements

- macOS (Apple Silicon M1/M2/M3) or Linux x86_64
- Docker Desktop (required for LinearDesign on macOS ARM)
- Python 3.11+
- 16GB+ RAM (MHCflurry models + proteome k-mer index)
- ~5GB disk (models + reference data)

## Python Dependencies

```bash
# Core
pip install pandas numpy biopython scipy pyyaml requests

# MHC Binding — MHCflurry 2.x (Apache 2.0)
pip install mhcflurry
mhcflurry-downloads fetch models_class1_presentation   # ~135MB download

# If mhcflurry-downloads not on PATH:
# /path/to/python/bin/mhcflurry-downloads fetch models_class1_presentation

# Immunogenicity — TensorFlow (for DeepImmuno-CNN)
# macOS: must use tensorflow 2.15.x (newer versions crash on macOS)
pip install "tensorflow==2.15.1"

# Protein sequence lookup
pip install pyensembl
python -c "
import pyensembl
ens = pyensembl.EnsemblRelease(release=110, species='human')
ens.download()
ens.index()
"
# Downloads ~1.5GB to ~/Library/Caches/pyensembl/

# Testing
pip install pytest
```

## External Tools

### DeepImmuno-CNN (immunogenicity prediction)
```bash
cd /path/to/project
git clone https://github.com/frankligy/DeepImmuno.git models/DeepImmuno

# Fix checkpoint files (dot-prefixed files cause issues with TF 2.15)
cd models/DeepImmuno/models/cnn_model_331_3_7
cp .data-00000-of-00001 model.data-00000-of-00001
cp .index model.index
cat > checkpoint <<'EOF'
model_checkpoint_path: "model"
all_model_checkpoint_paths: "model"
EOF
```

### LinearDesign (mRNA codon + structure optimization)
```bash
cd /path/to/project
git clone https://github.com/LinearDesignSoftware/LinearDesign.git

# On Linux x86_64: compile natively
cd LinearDesign && make

# On macOS ARM: the pre-compiled .so has ABI issues.
# LinearDesign runs via Docker instead (automatic in pipeline):
docker run --rm --platform linux/amd64 \
  -v "$(pwd)/LinearDesign":/build -w /build \
  gcc:11 make
```

**License**: Academic/research use free. Commercial requires separate license.

### Human Reference Proteome (foreignness scoring)
```bash
mkdir -p reference/proteome
curl -L -o reference/proteome/human_proteome.fasta.gz \
  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
gunzip -k reference/proteome/human_proteome.fasta.gz
# ~27MB uncompressed, 20,659 proteins
```

## Demo Data (TCGA-EE-A3J5 melanoma)

```bash
mkdir -p demo/tcga_skcm

# Somatic MAF (masked, open access)
curl -L -o demo/tcga_skcm/TCGA-EE-A3J5_somatic.maf.gz \
  "https://api.gdc.cancer.gov/data/3a7390dd-f069-4c6a-8802-a80f09662663"

# Gene expression (STAR counts with TPM)
curl -L -o demo/tcga_skcm/TCGA-EE-A3J5_expression.tsv \
  "https://api.gdc.cancer.gov/data/d076714e-3a5a-479e-81d5-2a377d3e37ea"
```

## Running the Pipeline

```bash
# Full end-to-end (VCF → mRNA)
python tests/run_e2e.py

# What it does:
# Step 0: MAF → peptides (Ensembl GRCh38) + MHCflurry binding predictions
# Module 8: DeepImmuno-CNN immunogenicity + foreignness (UniProt proteome)
# Module 9: Composite scoring + ranking + selection (top 20)
# Module 10: Polyepitope design (greedy TSP ordering, AAY linkers, MITD)
# Module 11: LinearDesign codon optimization + UTRs + poly-A

# Results saved to results/TCGA_EE_A3J5/
```

## Pipeline Tools Summary

| Step | Tool | License | Runs on |
|------|------|---------|---------|
| Protein sequences | pyensembl (Ensembl GRCh38 r110) | Apache 2.0 | Local |
| MHC-I binding | MHCflurry 2.2.0 | Apache 2.0 | Local |
| Immunogenicity | DeepImmuno-CNN | MIT | Local (TF) |
| Foreignness | k-mer vs UniProt proteome | - | Local |
| Ranking | Custom composite scoring | - | Local |
| Polyepitope | Greedy TSP + junctional QC | - | Local |
| Codon optimization | LinearDesign | Academic free | Docker (x86_64) |
| Structural scoring | AlphaFold2-Multimer (NVIDIA API) | API key required | Remote (async) |

## Known Limitations

- **CCF = 1.0 for all candidates**: PyClone-VI requires BAM files with VAF data. From MAF-only input, all mutations are assumed clonal.
- **No MHC Class II typing**: OptiType supports Class I only. Class II HLA typing requires HLA-HD (academic license).
- **Structural scoring not yet in e2e**: AlphaFold2-Multimer API is async (hours per job). Script exists but not wired into main pipeline run.
- **LinearDesign requires Docker on macOS ARM**: Pre-compiled x86_64 .so files have ABI mismatch with macOS ARM clang.

## Dockerization Notes (for server/cloud deployment)

The pipeline runs these tools:
1. **MHCflurry** — `pip install mhcflurry` + model download (~135MB)
2. **TensorFlow 2.15** — ~500MB
3. **DeepImmuno** — git clone + checkpoint files (~500KB)
4. **LinearDesign** — compiled C++ binary (Linux x86_64 only)
5. **pyensembl** — Ensembl data cache (~1.5GB)
6. **UniProt proteome** — ~27MB FASTA

Total Docker image size estimate: ~4GB base + ~3GB models/data = ~7GB

For AWS Batch / cloud: Nextflow profiles in `conf/aws.config` define queue mappings.
