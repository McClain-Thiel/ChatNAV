# Setup Log — Neoantigen Vaccine Pipeline

Everything needed to reproduce this environment from scratch.

## System Requirements

- macOS (Apple Silicon M1/M2/M3) or Linux x86_64
- Docker Desktop (required for LinearDesign on macOS ARM)
- Python 3.11+
- 16GB+ RAM (MHCflurry models + proteome k-mer index)
- ~5GB disk (models + reference data)
- GPU recommended for BigMHC-IM scoring (works on CPU but slower)

## Python Dependencies

```bash
# Core pipeline (exact versions pinned — see pyproject.toml)
pip install pandas==2.3.2 numpy==2.0.2 biopython==1.85 scipy==1.13.1 \
    pyyaml==6.0.2 requests==2.32.5 scikit-learn==1.6.1

# MHC Binding — MHCflurry 2.2.0 (Apache 2.0)
# IMPORTANT: pin both the package AND model release — different versions
# produce different percentile ranks
pip install mhcflurry==2.2.0
mhcflurry-downloads fetch models_class1_presentation   # model release 2.2.0, ~135MB

# Immunogenicity — BigMHC-IM 1.0.1 (primary, PyTorch-based)
pip install git+https://github.com/griffithlab/bigmhc.git@v1.0.1

# Protein sequence lookup — Ensembl release 110 (pinned)
pip install pyensembl
python -c "
import pyensembl
ens = pyensembl.EnsemblRelease(release=110, species='human')
ens.download()
ens.index()
"
# Downloads ~1.5GB to ~/Library/Caches/pyensembl/

# Benchmarking and ML reranker (optional)
pip install lightgbm

# Testing
pip install pytest
```

## External Tools

### BigMHC-IM (immunogenicity prediction — primary model)
```bash
pip install git+https://github.com/griffithlab/bigmhc.git
# Ensemble of 7 transformer models, runs on CPU or GPU
# No additional data downloads needed (models bundled in package)
```

### DeepImmuno-CNN (legacy, optional)
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
Note: DeepImmuno is no longer used in the default pipeline. BigMHC-IM replaced it.

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
# First run will auto-generate a k-mer cache (.kmers_k9.npy) for fast reloading
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
# Full end-to-end (MAF → mRNA)
python tests/run_e2e.py

# What it does:
# Step 0: MAF → peptides (Ensembl GRCh38) + MHCflurry binding predictions
# Module 8a: BigMHC-IM immunogenicity + foreignness (UniProt proteome)
# Module 8b: Structural scoring (TCR-facing position lookup)
# Module 9: Sequential filters + ranking + selection (top 20)
# Module 10: Polyepitope design (greedy TSP ordering, AAY linkers, MITD)
# Module 11: LinearDesign codon optimization + UTRs + poly-A

# Results saved to results/TCGA_EE_A3J5/
```

## Benchmarking

```bash
# Create patient splits (only needed once)
python benchmark/make_splits.py \
  --input benchmark/muller/Mutation_data_org.txt \
  --dev-frac 0.7 --seed 42 \
  --output benchmark/muller/splits.json

# Run benchmark on dev set
python benchmark/run_muller_benchmark.py \
  --split dev --bootstrap 1000 --output experiments/my_experiment/

# Compare against baseline
python benchmark/compare_runs.py \
  --baseline experiments/baseline_2026-04-07/ \
  --candidate experiments/my_experiment/
```

## Pipeline Tools Summary

| Step | Tool | License | Runs on |
|------|------|---------|---------|
| Protein sequences | pyensembl (Ensembl GRCh38 r110) | Apache 2.0 | Local |
| MHC-I binding | MHCflurry 2.2.0 | Apache 2.0 | Local |
| Immunogenicity | BigMHC-IM (ensemble of 7 models) | Academic | Local (PyTorch) |
| Foreignness | k-mer vs UniProt proteome | -- | Local |
| Ranking | LightGBM reranker or composite scoring | MIT | Local |
| Polyepitope | Greedy TSP + junctional QC | -- | Local |
| Codon optimization | LinearDesign | Academic free | Docker (x86_64) |
| Structural scoring | Position lookup / PANDORA / AlphaFold2 | Mixed | Local / Remote |

## Known Limitations

- **CCF = 1.0 for all candidates**: PyClone-VI requires BAM files with VAF data. From MAF-only input, all mutations are assumed clonal.
- **No MHC Class II typing**: OptiType supports Class I only. Class II HLA typing requires HLA-HD (academic license).
- **Structural scoring not yet in e2e**: AlphaFold2-Multimer API is async (hours per job). Script exists but not wired into main pipeline run.
- **LinearDesign requires Docker on macOS ARM**: Pre-compiled x86_64 .so files have ABI mismatch with macOS ARM clang.
- **Ensembl release pinned to 110**: Do not use `--release latest` — protein sequences change between releases.

## Dockerization Notes (for server/cloud deployment)

The pipeline runs these tools:
1. **MHCflurry** — `pip install mhcflurry` + model download (~135MB)
2. **BigMHC-IM** — PyTorch + 7 ensemble model weights (~200MB)
3. **LinearDesign** — compiled C++ binary (Linux x86_64 only)
4. **pyensembl** — Ensembl data cache (~1.5GB)
5. **UniProt proteome** — ~27MB FASTA + k-mer cache

Total Docker image size estimate: ~4GB base + ~2GB models/data = ~6GB

For AWS Batch / cloud: Nextflow profiles in `conf/aws.config` define queue mappings.
