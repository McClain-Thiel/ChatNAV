# ChatNAV — Personalized Neoantigen Vaccine Design Pipeline

An end-to-end computational pipeline for designing personalized mRNA neoantigen vaccines from somatic mutation data. Takes a patient's somatic variants (MAF or VEP-annotated VCF), HLA types, and gene expression as input and outputs a synthesis-ready mRNA vaccine sequence with confidence estimates.

---

## Background

### What is a neoantigen vaccine?

Cancer cells accumulate somatic mutations that produce altered proteins not found in normal tissue. When these mutant proteins are degraded inside the cell, short peptide fragments (8-11 amino acids for MHC Class I, 13-25 for Class II) are loaded onto MHC molecules and displayed on the cell surface. If a displayed peptide is sufficiently different from any normal human protein, circulating T cells can recognize it as foreign and mount an immune response. These tumor-specific peptides are called **neoantigens**.

A personalized neoantigen vaccine encodes a set of patient-specific neoantigens as a single mRNA molecule. When injected, the mRNA is translated by the patient's cells, processed through the same MHC presentation pathway, and displayed to T cells — training the immune system to recognize and kill tumor cells carrying those mutations.

### Clinical context

This approach is in late-stage clinical trials. BioNTech's **autogene cevumeran** (BNT122) and Moderna's **mRNA-4157** (V940) both encode ~20 neoantigens per patient as a polyepitope mRNA construct, individually designed and manufactured within weeks of tumor sequencing. In melanoma (BNT122 + anti-PD1, KEYNOTE-942) and pancreatic cancer (autogene cevumeran + atezolizumab), patients receiving personalized vaccines showed durable neoantigen-specific T cell responses and improved recurrence-free survival compared to checkpoint inhibitors alone.

### The computational problem

A typical tumor has hundreds to thousands of somatic mutations, but only ~1-5% produce peptides that are both presented by the patient's MHC alleles and actually immunogenic. The core challenge is **neoantigen selection**: choosing the best ~20 candidates from thousands.

This requires predicting multiple properties for each candidate:

1. **MHC binding**: Will the peptide bind the patient's specific HLA alleles tightly enough to be presented? (~2% of random peptides bind any given allele)
2. **Immunogenicity**: Even if presented, will the peptide trigger a T cell response? This depends on how foreign it looks vs the human self-proteome, the biochemistry of the peptide-MHC-TCR interaction, and whether the T cell repertoire contains receptors that recognize it.
3. **Expression**: Is the gene actually expressed in this patient's tumor? A mutation in a silenced gene produces no protein, so no peptide.
4. **Clonality**: Is the mutation present in all tumor cells (clonal) or only a subpopulation (subclonal)? Clonal mutations are better targets because every tumor cell displays them.

ChatNAV combines these signals using a LightGBM reranker trained on 131 patients with experimentally validated immunogenicity data, achieving Macro Recall@20 of 0.68 on held-out patients — meaning ~68% of a patient's validated immunogenic neoantigens appear in the pipeline's top 20 predictions.

### mRNA vaccine construct design

Once neoantigens are selected, the pipeline assembles them into a single mRNA construct following the architecture used by BNT122 and mRNA-4157:

```
5' Cap — 5' UTR — Signal Peptide — [Epitope-Linker]×20 — MITD — Stop — 3' UTR — Poly-A
```

- **Signal peptide** (21 aa): routes the polyepitope into the secretory pathway for MHC loading
- **AAY linkers** (MHC-I epitopes): cleavable by the proteasome, minimizing junctional neoepitopes
- **GPGPG linkers** (MHC-II epitopes): flexible spacers for Class II presentation
- **MITD domain** (58 aa): transmembrane domain from MHC-I that enhances antigen presentation
- **Codon optimization**: LinearDesign jointly optimizes codon usage (CAI) and mRNA secondary structure (MFE) for maximum translation efficiency and stability
- **N1-methylpseudouridine**: substituted at all uridine positions to reduce innate immune activation and increase mRNA half-life

---

## Pipeline: Engineering View

### What goes in

ChatNAV accepts **somatic variants** — the output of a somatic variant calling pipeline (e.g., Mutect2, Strelka2, nf-core/sarek). It does **not** do variant calling from raw sequencing data (BAM/FASTQ). The upstream workflow is:

```
Tumor BAM + Normal BAM
    → Somatic variant caller (Mutect2, Strelka2)
    → Somatic VCF (tumor-specific mutations only)
    → VEP annotation (adds gene, protein change, consequence)
    → ChatNAV starts here
```

**Required inputs (3 files + HLA alleles):**

| Input | Format | Source | What it provides |
|---|---|---|---|
| Somatic variants | MAF (`.maf`/`.maf.gz`) or VEP-annotated VCF (`.vcf`/`.vcf.gz`) | GDC, vcf2maf, sarek, Mutect2+VEP | Gene, protein change, variant class |
| Gene expression | TSV with gene_id + TPM columns | STAR, Salmon, Kallisto | Whether the mutated gene is expressed |
| HLA alleles | Comma-separated string | OptiType, HLA-HD, arcasHLA, clinical typing | Patient's MHC alleles for binding prediction |
| *(optional)* CCF | Column in MAF, or PyClone-VI output | PyClone-VI, ABSOLUTE, Sequenza | Clonality of each mutation |

**VCF note:** The `--vcf` flag accepts a single somatic VCF that already contains tumor-vs-normal calls (not two separate VCFs). If the VCF has multiple sample columns (e.g., TUMOR and NORMAL from Mutect2), use `--tumor-sample TUMOR` to specify which column to read allele depths from.

### What comes out

```bash
results/patient_001/
├── selected_neoantigens.tsv    # Top 20 epitopes with scores, HLA, flags
├── polyepitope.faa             # Polyepitope protein sequence
├── polyepitope_design.tsv      # Construct layout: epitope order, linkers, junction QC
├── mrna_sequence.fasta         # Synthesis-ready mRNA nucleotide sequence
├── synthesis_spec.json         # Full spec for synthesis facility (see below)
├── summary_card.html           # One-page visual report
├── confidence.json             # Vaccine potency estimate
├── binding_predictions.tsv     # All MHC binding predictions (intermediate)
├── candidates_meta.tsv         # All candidate peptides with metadata (intermediate)
├── immunogenicity_scores.tsv   # BigMHC-IM + foreignness scores (intermediate)
└── structural_scores.tsv       # Structural tier scores (intermediate)
```

**Key deliverables:**

- **`synthesis_spec.json`** — everything a synthesis facility needs: full mRNA sequence, CDS boundaries, UTR sequences, pseudouridine positions, poly-A length, cap analog, GC content, codon optimization stats (MFE, CAI), and CDS quality scans (internal poly-A signals, cryptic splice sites)
- **`selected_neoantigens.tsv`** — the 20 selected epitopes with: gene, mutation, peptide sequence, HLA allele, MHC class (I/II), binding rank, immunogenicity score, foreignness, agretopicity, validation flags (e.g., `high_immunogenicity; strong_differential_binding; very_stable_pMHC`)
- **`confidence.json`** — empirical potency estimate based on mutation count vs benchmark performance: confidence tier (HIGH/MODERATE/LOW), expected Recall@20, and caveats

### How it runs

```bash
# From MAF
python bin/run_pipeline.py \
    --maf patient_somatic.maf.gz \
    --expression patient_expression.tsv \
    --hla-alleles 'HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*44:03,HLA-C*07:02,HLA-C*05:01' \
    --patient-id PT-001 \
    --output results/PT-001/

# From VEP-annotated somatic VCF
python bin/run_pipeline.py \
    --vcf patient_somatic.vep.vcf.gz \
    --tumor-sample TUMOR \
    --expression patient_expression.tsv \
    --hla-alleles 'HLA-A*02:01,HLA-B*07:02,HLA-C*07:02' \
    --patient-id PT-002 \
    --output results/PT-002/

# Validate inputs without running (checks files, HLA format, tool availability)
python bin/run_pipeline.py --maf patient.maf --expression expr.tsv \
    --hla-alleles 'HLA-A*02:01' --output /tmp/test --dry-run
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `--profile` | `high_tmb` | Scoring profile: `high_tmb` (melanoma, NSCLC), `low_tmb` (PDAC, GBM), `research` (no filters) |
| `--top-n` | 20 | Number of epitopes to select |
| `--max-mutations` | all | Limit mutations processed (for testing) |
| `--skip-mrna` | false | Skip LinearDesign step (useful if Docker unavailable) |
| `--dry-run` | false | Validate inputs only |

### Pipeline stages and timing

| Stage | Script | What it does | Time (50 mutations) |
|---|---|---|---|
| 0a | `vcf_to_maf.py` | VCF → MAF conversion (if VCF input) | <1s |
| 0b | `maf_to_pipeline_input.py` | MAF → peptide windows + MHCflurry binding + MHCnuggets Class II | ~30s |
| 8a | `score_immunogenicity.py` | BigMHC-IM immunogenicity + k-mer foreignness + agretopicity | ~40s (CPU), ~5s (GPU) |
| 8b | `structural_scoring.py` | TCR-facing position lookup | <1s |
| 9 | `rank_and_select.py` | Hard filters → LightGBM reranker → class balance → HLA diversity | <1s |
| 10 | `design_polyepitope.py` | Greedy TSP ordering + linkers + junction QC | ~3s |
| 11 | `design_mrna.py` | LinearDesign codon optimization + UTRs + poly-A + CDS QC scans | ~5s |
| -- | `generate_summary_card.py` | HTML summary report | <1s |
| | | **Total** | **~75s** (CPU) |

BigMHC-IM is the bottleneck. On a GPU (L40S), scoring drops from ~40s to ~5s for 50 mutations.

---

## Pipeline: Biological View

### Stage 1: Candidate generation

For each somatic mutation (missense, frameshift, in-frame indel), the pipeline:

1. **Looks up the protein sequence** from Ensembl GRCh38 (release 110) using the gene symbol and protein coordinates from the MAF/VCF
2. **Generates all peptide windows** that contain the mutation site: 8-11mers for MHC Class I, 13-17mers for MHC Class II
3. **Predicts MHC binding** using MHCflurry 2.2 (Class I, Apache 2.0) and MHCnuggets (Class II, MIT) for each peptide-HLA pair. Also predicts wildtype peptide binding at the same positions for agretopicity calculation

A typical patient with 500 mutations produces ~15,000 peptide-HLA pairs.

### Stage 2: Immunogenicity scoring

Each candidate peptide is scored on multiple orthogonal signals:

- **BigMHC-IM** (importance rank: #1): ensemble of 7 transformer models trained on mass spectrometry immunopeptidomics data. Predicts whether a peptide-MHC complex will be immunogenic (i.e., trigger a T cell response), beyond just binding.
- **Expression** (#2): `log10(TPM + 1)` normalized to [0, 1]. A mutation in a silenced gene cannot produce a neoantigen.
- **HydroCore** (#3, IMPROVE feature): mean Kyte-Doolittle hydrophobicity of the peptide's binding core positions (P4-P7 for 9-mers). Hydrophobic residues at TCR-contact positions correlate with immunogenicity.
- **Agretopicity** (#4): `log2(wt_rank / mut_rank)`. High agretopicity means the mutation creates a new MHC binder that the wildtype doesn't have — the immune system hasn't seen it before, so tolerance is unlikely.
- **CCF** (#5): cancer cell fraction. Clonal mutations (CCF ~1.0) are present in every tumor cell; subclonal mutations may be lost to immune selection.
- **Binding strength** (#6): MHCflurry presentation percentile rank. Strong binders (< 2%) are more likely to be displayed.
- **Foreignness**: k-mer overlap against 10.4M 9-mers from the UniProt human proteome. Peptides dissimilar to all human proteins are more likely to be recognized as non-self.
- **Differential agretopicity**: binary flag — mutant binds well (< 0.5%) AND wildtype doesn't (> 2%). Strong differential presentation.
- **PropHydroAro** (IMPROVE): proportion of hydrophobic + aromatic residues (FYWILVM) in the binding core.

### Stage 3: Ranking and selection

1. **Hard filters** remove candidates that fail biological prerequisites:
   - MHC binding rank > 2% (not presented)
   - Gene TPM < 1 (not expressed)
   - CCF < 0.3 (too subclonal)
   - Exact self-match in human proteome (tolerance risk)

2. **LightGBM reranker** scores surviving candidates using 9 features (trained on 91 patients with experimentally validated immunogenicity). Candidates are ranked by predicted probability of immunogenicity.

3. **Class balance**: 70-85% of slots go to MHC Class I epitopes (CD8+ T cell targets), 15-30% to Class II (CD4+ T cell help). Both are needed for durable responses — Class II epitopes provide helper signals that sustain the CD8+ response.

4. **HLA diversity**: ensures the top 20 covers multiple HLA alleles, not just the easiest-to-predict one.

5. **Frameshift and shared neoantigen bonuses**: frameshifts produce entirely novel protein sequences (maximum foreignness), and shared neoantigens (e.g., KRAS G12V) have confirmed immunogenicity across patients.

### Stage 4: Construct design

The 20 selected epitopes are assembled into a polyepitope protein:

- **Ordering**: greedy TSP minimizes junctional neoepitopes — peptides are arranged so that the junction between adjacent epitopes + linker does not accidentally create a strong MHC binder (which would compete with the intended epitopes for presentation)
- **Junction QC**: every 8-11mer spanning a linker junction is checked for MHC binding. Junctions with binding rank < 2% are flagged.
- **mRNA design**: the protein is reverse-translated with LinearDesign (joint codon usage + mRNA structure optimization), wrapped in UTRs and poly-A tail, and scanned for internal poly-A signals and cryptic splice sites that could cause premature termination

---

## Benchmark Results

### Muller/Gfeller Harmonized Dataset (Immunity 2023)

131 patients with experimentally validated CD8+ T cell responses to individual neoantigens. Split into dev (91) and validation (39) sets. All results use macro Recall@20 with 1000-bootstrap patient-resampled confidence intervals.

**Pipeline-compatible reranker (9 features, used in production):**

| Metric | Hand-tuned Baseline | LightGBM Reranker |
|---|---|---|
| Dev Macro Recall@20 | 0.414 [0.310, 0.521] | 0.635 [0.547, 0.724] |
| **Val Macro Recall@20** | 0.517 [0.353, 0.693] | **0.684 [0.547, 0.814]** |
| Val AUC | — | 0.750 |

**Feature importance (9-feature production model):**

| Feature | LightGBM Gain | What it captures |
|---|---|---|
| HydroCore | 409 | Binding core hydrophobicity → TCR recognition |
| BigMHC-IM | 386 | Learned immunogenicity from mass spec data |
| Expression | 367 | Gene must be expressed to produce peptide |
| CCF | 321 | Clonal mutations target all tumor cells |
| Agretopicity | 287 | Differential binding vs wildtype |
| Binding strength | 240 | Presentation probability |
| PropHydroAro | 58 | Aromatic/hydrophobic content in core |
| Foreignness | 38 | Dissimilarity to human self-proteome |
| Diff agretopicity | 1 | Binary strong-differential flag |

See `experiments/LOG.md` for the full experiment history (22 experiments across 8 phases).

---

## Quick Start

```bash
# 1. Install dependencies (exact versions pinned in pyproject.toml)
pip install pandas==2.3.2 numpy==2.0.2 biopython==1.85 scipy==1.13.1 \
    pyyaml==6.0.2 requests==2.32.5 scikit-learn==1.6.1 pyensembl==2.3.13 \
    mhcflurry==2.2.0 lightgbm==4.6.0
pip install git+https://github.com/griffithlab/bigmhc.git

# Download MHCflurry models (~135MB)
mhcflurry-downloads fetch models_class1_presentation

# Download Ensembl reference (~1.5GB, cached)
python -c "
import pyensembl
ens = pyensembl.EnsemblRelease(110, species='human')
ens.download(); ens.index()
"

# 2. Initialize submodules
git submodule update --init --recursive

# 3. Download reference proteome
mkdir -p reference/proteome
curl -L -o reference/proteome/human_proteome.fasta.gz \
  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
gunzip -k reference/proteome/human_proteome.fasta.gz

# 4. Build LinearDesign (requires Docker on macOS ARM, native on Linux)
cd LinearDesign && make; cd ..

# 5. Run on demo patient (TCGA melanoma)
python bin/run_pipeline.py \
    --maf demo/tcga_skcm/TCGA-EE-A3J5_somatic.maf.gz \
    --expression demo/tcga_skcm/TCGA-EE-A3J5_expression.tsv \
    --hla-alleles 'HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*44:03,HLA-C*07:02,HLA-C*05:01' \
    --patient-id TCGA-EE-A3J5 \
    --output results/TCGA-EE-A3J5/
```

---

## Project Structure

```
ChatNAV/
├── bin/                            # Pipeline scripts
│   ├── run_pipeline.py                 # Unified entry point (MAF/VCF → vaccine)
│   ├── vcf_to_maf.py                  # VEP-annotated VCF → MAF conversion
│   ├── maf_to_pipeline_input.py       # MAF → peptides + MHCflurry/MHCnuggets
│   ├── score_immunogenicity.py        # BigMHC-IM + foreignness + agretopicity
│   ├── structural_scoring.py          # TCR-facing position lookup
│   ├── rank_and_select.py             # Filters → LightGBM reranker → selection
│   ├── design_polyepitope.py          # Polyepitope assembly + junction QC
│   ├── design_mrna.py                 # LinearDesign + UTRs + CDS QC scans
│   └── generate_summary_card.py       # HTML vaccine report
├── models/
│   ├── reranker/                      # Production LightGBM model (9 features)
│   │   ├── lgbm_pipeline_v1.txt       # Frozen model weights
│   │   └── pipeline_features.json     # Feature names + mapping from pipeline columns
│   ├── DeepImmuno/                    # Legacy immunogenicity model (not used)
│   └── CHECKSUMS.txt                  # SHA256 of all model files
├── benchmark/                      # Benchmark data + evaluation
│   ├── muller/                        # Muller/Gfeller 131-patient dataset
│   ├── gartner/                       # Internal validation (touch-once)
│   ├── tesla/                         # Final holdout (release-only)
│   ├── run_muller_benchmark.py        # Benchmark runner with bootstrap CIs
│   └── compare_runs.py               # CI overlap comparison
├── experiments/                    # Experiment log + artifacts
│   └── LOG.md                         # 22 experiments, append-only
├── conf/scoring_weights.yaml       # Filter thresholds + profile configs
├── reference/                      # Proteome, UTRs, signal peptides
├── docker/                         # 14 Dockerfiles (one per stage)
├── tests/                          # 55 tests (pytest + e2e)
├── INFRASTRUCTURE.md               # Deployment architecture + service roadmap
├── SETUP.md                        # Detailed installation guide
└── plan.md                         # System design document
```

## Infrastructure

See [INFRASTRUCTURE.md](INFRASTRUCTURE.md) for the full deployment architecture, AWS resources, and the roadmap for a FastAPI service + web frontend.

- **GPU instance**: g6e.4xlarge (NVIDIA L40S 46GB) for BigMHC scoring
- **S3**: `s3://chatnav-pipeline-data/` (~620GB cached reference data)
- **14 Docker images**: one per pipeline stage, in `docker/`

## Known Limitations

- **No variant calling**: pipeline requires pre-called somatic variants (MAF or VEP-annotated VCF). Use Mutect2/Strelka2/sarek upstream.
- **CCF defaults to 1.0**: from MAF-only input without PyClone-VI, all mutations are assumed clonal. Real CCF requires BAM files with VAF data.
- **No MHC Class II typing**: OptiType supports Class I only. Class II alleles must be provided manually (from HLA-HD or clinical typing).
- **LinearDesign requires Docker on macOS ARM**: pre-compiled binary is x86-only. Use `--skip-mrna` to bypass, or run on Linux.
- **Structural scoring demoted**: Tier 1 position lookup (AUC 0.489) was found to add no value in benchmark ablation (EXP-040). Available in `research` profile only.

## References

- Muller et al. (2023). Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction. *Immunity*.
- Albert et al. (2023). BigMHC: accurate prediction of MHC-I antigen presentation and immunogenicity. *Nature Machine Intelligence*.
- O'Donnell et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides. *Cell Systems*.
- Zhang et al. (2023). LinearDesign: Efficient algorithms for optimized mRNA sequence design. *Nature*.
- Gartner et al. (2021). A machine learning model for ranking candidate neoantigens. *Nature Cancer*.
- Wells et al. (2020). Key Parameters of Tumor Epitope Immunogenicity. *Cell*.
- Sahin et al. (2017). Personalized RNA mutanome vaccines mobilize poly-specific therapeutic immunity against cancer. *Nature*.
- Rojas et al. (2023). Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer. *Nature*.

## License

Pipeline code: MIT. External tools (BigMHC, LinearDesign, MHCflurry) have their own licenses — see SETUP.md.
