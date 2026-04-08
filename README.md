# ChatNAV — Personalized Neoantigen Vaccine Design Pipeline

An end-to-end computational pipeline for designing personalized mRNA neoantigen vaccines from somatic mutation data. Takes a patient's somatic mutations (MAF/VCF), HLA types, and gene expression as input and outputs a synthesis-ready mRNA vaccine sequence.

## Background

Cancer cells accumulate somatic mutations that produce altered proteins not found in normal tissue. Short peptide fragments (8-11 amino acids) from these mutant proteins can be presented on the cell surface by MHC molecules, where they are visible to the immune system as **neoantigens**. A personalized neoantigen vaccine encodes a set of patient-specific neoantigens as a single mRNA construct, training the patient's T cells to recognize and kill tumor cells displaying those peptides.

This approach is already in late-stage clinical trials. BioNTech's **autogene cevumeran** (BNT122) and Moderna's **mRNA-4157** (V940) both use individualized mRNA vaccines encoding ~20 neoantigens per patient, with each vaccine designed and manufactured within weeks of tumor sequencing. Early results in melanoma and pancreatic cancer show durable T cell responses and improved recurrence-free survival when combined with checkpoint inhibitors.

The computational challenge is **neoantigen selection**: a typical tumor has hundreds to thousands of somatic mutations, but only a small fraction (~1-5%) produce peptides that are (a) presented by the patient's MHC alleles and (b) actually immunogenic. Selecting the right 20 neoantigens from thousands of candidates is the core prediction problem this pipeline addresses.

ChatNAV implements the full pipeline from raw mutation calls to a synthesis-ready mRNA sequence, with a focus on principled neoantigen ranking. The current best model (a LightGBM reranker trained on 13 immunogenicity-related features) achieves Macro Recall@20 of 0.60 on the Muller/Gfeller benchmark dataset, meaning that on average 60% of a patient's validated immunogenic neoantigens appear in the top 20 predictions.

## Pipeline Overview

```
Somatic MAF/VCF + Expression + HLA alleles
    │
    ▼
Step 0: Peptide Generation + MHC Binding
    Ensembl GRCh38 protein lookup → 8-11mer sliding windows
    MHCflurry 2.2 presentation predictions (mutant + wildtype)
    │
    ▼
Module 8a: Scoring Signals
    BigMHC-IM immunogenicity score
    k-mer foreignness vs human proteome
    Agretopicity (wt vs mut binding ratio)
    │
    ▼
Module 8b: Structural Scoring (Tiered)
    Tier 1: TCR-facing position lookup — instant
    Tier 2: PANDORA homology modeling + SASA — minutes, CPU
    Tier 3: AlphaFold2-Multimer (local NIM) — hours, GPU
    │
    ▼
Module 9: Sequential Filter Pipeline
    1. FILTER: binding rank < 2%
    2. FILTER: expression TPM > 1
    3. FILTER: CCF > 0.3 (clonal)
    4. FILTER: not self-similar
    5. RANK: LightGBM reranker (13 features) or binding strength
    6. BONUS: frameshifts + shared neoantigens promoted
    7. SELECT: top 20 with HLA diversity
    │
    ▼
Module 10: Polyepitope Design
    Greedy TSP ordering (minimize junctional epitopes)
    Signal peptide + [Epitope1—AAY—Epitope2—AAY—...] + MITD domain
    │
    ▼
Module 11: mRNA Design
    LinearDesign: joint codon + mRNA structure optimization
    5'cap — 5'UTR — CDS — 3'UTR — poly-A → synthesis_spec.json
    │
    ▼
Output: synthesis-ready mRNA sequence + specification
```

## Final Output

The pipeline produces a single mRNA construct encoding all selected neoantigens:

```
5' Cap — 5' UTR (HBB, 45nt) — CDS — 3' UTR (AES/mtRNR1, 119nt) — Poly-A (120nt)
                                 │
                    Signal Peptide (21 aa)
                    Epitope 1 — AAY linker
                    Epitope 2 — AAY linker
                    ...
                    Epitope 20 — AAAAA
                    MITD domain (58 aa)
                    Stop codon
```

- **Codon optimized** by LinearDesign (joint CAI + MFE optimization)
- **N1-methylpseudouridine** substitution at all uridine positions
- **CleanCap-AG** co-transcriptional capping
- Output files: `mrna_sequence.fasta` + `synthesis_spec.json`

This is the format used by BioNTech (BNT122) and Moderna (mRNA-4157) for their individualized neoantigen vaccines.

## Tools

| Step | Tool | What it does | License |
|------|------|-------------|---------|
| Protein sequences | pyensembl (Ensembl GRCh38) | Real protein sequences from gene + mutation annotation | Apache 2.0 |
| MHC-I binding | MHCflurry 2.2 | Binding affinity + antigen processing + presentation | Apache 2.0 |
| Immunogenicity | BigMHC-IM | Transformer-based immunogenicity prediction | Academic |
| Foreignness | k-mer vs UniProt proteome | Self-similarity scoring (10.4M 9-mers) | -- |
| Structural | Position lookup / PANDORA / AlphaFold2 | TCR exposure scoring (tiered) | Mixed |
| Codon optimization | LinearDesign | Joint codon + mRNA structure optimization | Academic |

## Quick Start

```bash
# 1. Install dependencies
pip install mhcflurry pyensembl pandas numpy biopython pyyaml requests
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

# 4. Build LinearDesign (requires Docker on macOS ARM)
cd LinearDesign && make
cd ..

# 5. Download demo data (TCGA melanoma patient)
mkdir -p demo/tcga_skcm
curl -L -o demo/tcga_skcm/TCGA-EE-A3J5_somatic.maf.gz \
  "https://api.gdc.cancer.gov/data/3a7390dd-f069-4c6a-8802-a80f09662663"
curl -L -o demo/tcga_skcm/TCGA-EE-A3J5_expression.tsv \
  "https://api.gdc.cancer.gov/data/d076714e-3a5a-479e-81d5-2a377d3e37ea"

# 6. Run the pipeline
python bin/run_pipeline.py \
    --maf demo/tcga_skcm/TCGA-EE-A3J5_somatic.maf.gz \
    --expression demo/tcga_skcm/TCGA-EE-A3J5_expression.tsv \
    --hla-alleles 'HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*44:03,HLA-C*07:02,HLA-C*05:01' \
    --patient-id TCGA-EE-A3J5 \
    --output results/TCGA-EE-A3J5/
```

### Pipeline outputs

| File | Description |
|---|---|
| `selected_neoantigens.tsv` | Ranked epitope list with all scores and validation flags |
| `polyepitope.faa` | Polyepitope amino acid sequence (signal peptide + epitopes + linkers + MITD) |
| `polyepitope_design.tsv` | Construct layout with linker types and junction QC scores |
| `mrna_sequence.fasta` | Synthesis-ready mRNA (codon-optimized, UTRs, poly-A) |
| `synthesis_spec.json` | Full specification: sequence, pseudouridine positions, QC scans |
| `summary_card.html` | One-page visual report for review |
| `confidence.json` | Expected vaccine potency estimate based on benchmark data |

## How Each Step Works

### Step 0: Peptide Generation + MHC Binding

Parses somatic mutations from TCGA masked MAF format. For each missense, frameshift, and in-frame indel:
1. Looks up the real protein sequence from Ensembl GRCh38 (release 110)
2. Generates all 8-11mer peptide windows containing the mutation site
3. Runs MHCflurry Class1PresentationPredictor for binding + presentation scores
4. Runs MHCflurry on wildtype peptides at the same positions for agretopicity

### Module 8a: Scoring Signals

Three independent signals computed for each peptide-HLA pair:

- **BigMHC-IM**: Transformer-based immunogenicity prediction (ensemble of 7 models). Used as a validation flag on top candidates, not as a primary ranking signal.
- **Foreignness**: k-mer overlap against the full UniProt human reference proteome. Score of 0 = exact self-match, 1 = maximally foreign.
- **Agretopicity**: `log2(wt_binding_rank / mut_binding_rank)`. Positive = mutation creates a new binding event. Negative = wildtype already binds (tolerance risk).

### Module 8b: Structural Scoring

Tiered approach to score whether the mutated residue is TCR-exposed:

- **Tier 1** (instant): Position lookup table from crystal structure surveys. P4-P7 are TCR-facing for 9-mers.
- **Tier 2** (minutes, CPU): PANDORA homology modeling + BioPython SASA computation. Requires MODELLER (free academic).
- **Tier 3** (hours, GPU): AlphaFold2-Multimer via local NVIDIA NIM container. Full pMHC structure prediction with real SASA scoring. Deployed on g6e.4xlarge (L40S 46GB).

### Module 9: Sequential Filter Pipeline

Sequential hard filters reduce the candidate pool, then candidates are ranked and selected:

1. **Binding filter**: MHCflurry presentation percentile < 2% (keeps strong binders)
2. **Expression filter**: Gene TPM > 1 (gene must be expressed in the tumor)
3. **CCF filter**: Cancer cell fraction > 0.3 (mutation must be clonal/subclonal)
4. **Self-match filter**: Exclude peptides identical to human proteome
5. **Rank**: By expression-weighted binding strength (default) or LightGBM reranker using 13 features including BigMHC-IM immunogenicity, IMPROVE hydrophobicity features, agretopicity, driver gene status, and TCR-facing position
6. **Bonus promotion**: Frameshifts and shared neoantigens get rank boosts
7. **HLA diversity**: Ensure coverage across patient's HLA alleles
8. **Select top 20**

### Module 10: Polyepitope Design

Assembles selected neoantigens into an optimized polyepitope construct:
- Greedy TSP ordering to minimize junctional neoepitope formation
- AAY linkers (MHC-I) / GPGPG linkers (MHC-II)
- Signal peptide for secretory pathway + MITD domain for MHC-I loading
- Junctional QC flags high-risk junctions

### Module 11: mRNA Design

Converts the polyepitope protein into a synthesis-ready mRNA using **LinearDesign** for joint optimization of codon usage (CAI) and mRNA secondary structure (MFE). Adds 5' UTR, 3' UTR, poly-A tail, and outputs a synthesis specification with pseudouridine positions.

## Example Output (TCGA-EE-A3J5 Melanoma)

```
Patient: TCGA-EE-A3J5 (TCGA SKCM melanoma)
Variants: 2,584 somatic → 1,890 peptide windows → 260 strong binders
  → 139 pass all filters → 20 selected
Polyepitope: 332 amino acids
mRNA: 1,253 nt | GC: 57.2% | MFE: -679.9 kcal/mol | CAI: 0.900
```

## Benchmark Results

### Muller/Gfeller Harmonized Dataset (Immunity 2023)

Primary benchmark: 131 patients split into dev (91) and validation (39) sets. Evaluated on 8,891 screened mutations (146 immunogenic) in the dev set. Bootstrap CIs by resampling patients (not peptides).

**Current best (EXP-013: LightGBM reranker, dev set):**

| Metric | Hand-tuned Composite | LightGBM Reranker |
|---|---|---|
| Macro Recall@20 | 0.414 [0.310, 0.521] | **0.599 [0.494, 0.694]** |
| Macro Recall@50 | 0.767 [0.671, 0.852] | **0.886 [0.822, 0.943]** |

The reranker uses 13 features with leave-one-patient-out cross-validation.

**Feature importance (LightGBM gain):**

| Feature | Importance | AUC (individual) |
|---|---|---|
| BigMHC-IM | 390 | 0.562 |
| Expression (log TPM) | 388 | 0.737 |
| HydroCore (IMPROVE) | 313 | 0.540 |
| Agretopicity | 299 | 0.594 |
| CCF | 274 | 0.556 |
| PRIME binding rank | 224 | 0.735 |
| Structural tier 1 | 89 | 0.489 |
| Driver mutation (Intogen) | 62 | 0.531 |
| PropHydroAro (IMPROVE) | 58 | 0.547 |
| Foreignness (k-mer) | 27 | 0.478 |

See `experiments/LOG.md` for the full experiment history and `experiments/baseline_2026-04-07/` for baseline artifacts.

## Project Structure

```
ChatNAV/
├── bin/                            # Pipeline scripts (one per stage)
│   ├── maf_to_pipeline_input.py        # Step 0: MAF → peptides + MHCflurry
│   ├── score_immunogenicity.py         # Module 8a: BigMHC-IM + foreignness
│   ├── structural_scoring.py           # Module 8b: Tiered structural scoring
│   ├── rank_and_select.py              # Module 9: Sequential filter ranking
│   ├── design_polyepitope.py           # Module 10: Polyepitope assembly
│   └── design_mrna.py                  # Module 11: mRNA + LinearDesign
├── benchmark/                      # Benchmark data + evaluation scripts
│   ├── muller/                         # Muller/Gfeller 131-patient dataset
│   │   ├── splits.json                 # Frozen dev/val split (91/39, seed=42)
│   │   └── Mutation_data_org.txt       # Raw mutation data
│   ├── gartner/                        # Gartner NCI dataset (internal validation)
│   ├── tesla/                          # TESLA dataset (final holdout)
│   ├── make_splits.py                  # Create patient-stratified splits
│   ├── run_muller_benchmark.py         # Full benchmark with bootstrap CIs
│   └── compare_runs.py                 # Compare experiments via CI overlap
├── conf/                           # Configuration
│   ├── scoring_weights.yaml            # Filter thresholds + profiles
│   └── exp*.yaml                       # Per-experiment config overrides
├── experiments/                    # Experiment outputs + log (append-only)
│   ├── LOG.md                          # Experiment log (source of truth)
│   ├── baseline_2026-04-07/            # Frozen baseline artifacts
│   └── exp*/                           # Per-experiment artifact dirs
├── reference/                      # Reference data
│   ├── proteome/                       # UniProt human proteome + k-mer cache
│   ├── shared_neoantigens/
│   └── utr_library/
├── LinearDesign/                   # Submodule: codon optimization
├── models/DeepImmuno/              # Submodule: legacy immunogenicity model
├── docker/                         # Dockerfiles per module
├── modules/                        # Nextflow DSL2 modules
├── subworkflows/                   # Nextflow subworkflows
├── tests/                          # pytest suite + e2e smoke test
├── main.nf                         # Nextflow entry point
├── nextflow.config
├── SETUP.md                        # Installation guide
└── plan.md                         # System design document
```

## Infrastructure

Current deployment is a Nextflow pipeline running locally or on AWS Batch. See [INFRASTRUCTURE.md](INFRASTRUCTURE.md) for the full architecture, AWS resources, and the plan for wrapping this as a FastAPI service with a web frontend.

- **S3 bucket**: `s3://chatnav-pipeline-data/` (~620GB cached data)
- **GPU instance**: g6e.4xlarge (NVIDIA L40S 46GB) for BigMHC scoring and benchmarking
- **14 Docker images**: one per pipeline stage, in `docker/`
- New instances: `bash scripts/restore_from_s3.sh` pulls all cached data in ~5 minutes

## Known Limitations

- **CCF from MAF only**: Cancer cell fraction requires PyClone-VI with BAM files. From MAF-only input, CCF defaults to 1.0.
- **No MHC Class II typing**: OptiType supports Class I only.
- **MHCflurry vs NetMHCpan**: Uses MHCflurry (Apache 2.0). The 2% threshold is calibrated against MHCflurry's percentile distribution.
- **AlphaFold2 MSA search**: Takes ~45 min per prediction due to jackhmmer against large databases. Use `predict-structure-from-msa` endpoint with pre-computed MSAs for faster inference.

## References

- Albert et al. (2023). BigMHC: accurate prediction of MHC-I antigen presentation and immunogenicity. *Nature Machine Intelligence*.
- Muller et al. (2023). Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction. *Immunity*.
- O'Donnell et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides. *Cell Systems*.
- Zhang et al. (2023). LinearDesign: Efficient algorithms for optimized mRNA sequence design. *Nature*.
- Gartner et al. (2021). A machine learning model for ranking candidate neoantigens. *Nature Cancer*.
- Wells et al. (2020). Key Parameters of Tumor Epitope Immunogenicity. *Cell*.

## License

Pipeline code: MIT. External tools (BigMHC, LinearDesign, MHCflurry) have their own licenses -- see SETUP.md.
