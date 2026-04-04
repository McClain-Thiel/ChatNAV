# ChatNAV — Personalized Neoantigen Vaccine Design Pipeline

An end-to-end computational pipeline for designing personalized mRNA neoantigen vaccines from somatic mutation data. Takes a patient's somatic mutations (MAF/VCF), HLA types, and gene expression as input and outputs a synthesis-ready mRNA vaccine sequence.

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
    5. RANK: by binding strength
    6. BONUS: frameshifts + shared neoantigens promoted
    7. FLAG: immunogenicity, structural, agretopicity as validation
    8. SELECT: top 20 with HLA diversity
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
python tests/run_e2e.py
```

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

Sequential hard filters reduce the candidate pool, then ranking by binding strength:

1. **Binding filter**: MHCflurry presentation percentile < 2% (keeps strong binders)
2. **Expression filter**: Gene TPM > 1 (gene must be expressed in the tumor)
3. **CCF filter**: Cancer cell fraction > 0.3 (mutation must be clonal/subclonal)
4. **Self-match filter**: Exclude peptides identical to human proteome
5. **Rank by binding**: Strongest binders first within the filtered set
6. **Bonus promotion**: Frameshifts and shared neoantigens get rank boosts
7. **Validation flags**: Each candidate tagged with immunogenicity, structural, agretopicity assessments
8. **HLA diversity**: Ensure coverage across patient's HLA alleles
9. **Select top 20**

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

Primary benchmark: 131 patients, 13,483 screened mutations (213 immunogenic), all signals real (expression, CCF, agretopicity, wildtype sequences).

**Sequential filter funnel:**
```
13,483 mutations (213 immunogenic, 1.6%)
  ↓ Binding < 2%
9,369 (196 pos, 92% recall)
  ↓ Expression > 1 TPM
6,513 (184 pos, 86% recall)
  ↓ CCF > 0.3
~6,400 (180 pos, 85% recall)
  ↓ Rank by binding strength → top 20 per patient
```

**Per-patient results:**
| Metric | Value |
|---|---|
| Macro Recall@20 | **0.584** |
| Macro Recall@50 | 0.649 |
| Patients with R@20 > 0 | 72/97 (74%) |

**Individual signal AUCs (mutation-level):**
| Signal | AUC | Role in pipeline |
|---|---|---|
| PRIME binding rank | **0.746** | Primary ranking signal |
| Expression (TPM) | **0.739** | Hard filter |
| Agretopicity | 0.590 | Validation flag |
| BigMHC-IM | 0.586 | Validation flag |
| CCF | 0.535 | Hard filter |
| Structural tier 1 | 0.500 | Validation flag |
| Foreignness | 0.478 | Hard filter |

**Key finding**: Sequential hard filters + rank by binding strength (Recall@20 = 0.584) outperforms a weighted composite of all signals (Recall@20 = 0.372) by 57%. Binding rank is the strongest discriminator after filtering; immunogenicity and structural scores are most valuable as validation flags on top candidates.

## Project Structure

```
ChatNAV/
├── bin/                            # Pipeline scripts
│   ├── maf_to_pipeline_input.py        # Step 0: MAF → peptides + MHCflurry
│   ├── score_immunogenicity.py         # Module 8a: BigMHC-IM + foreignness
│   ├── structural_scoring.py           # Module 8b: Tiered structural scoring
│   ├── rank_and_select.py              # Module 9: Sequential filter ranking
│   ├── design_polyepitope.py           # Module 10: Polyepitope assembly
│   └── design_mrna.py                  # Module 11: mRNA + LinearDesign
├── conf/                           # Configuration
│   └── scoring_weights.yaml            # Filter thresholds + profiles
├── reference/                      # Small reference data
│   ├── signal_peptides.fasta
│   ├── shared_neoantigens/
│   └── utr_library/
├── benchmark/                      # Benchmark data + scripts
│   ├── muller/                         # Muller/Gfeller 131-patient dataset
│   └── gartner/                        # Gartner NCI 61-patient dataset
├── LinearDesign/                   # Submodule: codon optimization
├── models/DeepImmuno/              # Submodule: legacy model
├── docker/                         # Dockerfiles per module
├── modules/                        # Nextflow DSL2 modules
├── subworkflows/                   # Nextflow subworkflows
├── scripts/
│   └── restore_from_s3.sh             # Restore cached data from S3
├── tests/
│   └── run_e2e.py                      # Full e2e test
├── main.nf                         # Nextflow entry point
├── nextflow.config
├── SETUP.md                        # Detailed installation log
└── plan.md                         # System design document
```

## Infrastructure

- **S3 bucket**: `s3://chatnav-pipeline-data/` (~620GB cached data)
- **GPU instance**: g6e.4xlarge (NVIDIA L40S 46GB), AlphaFold2-Multimer NIM deployed
- **PANDORA**: Docker container with MODELLER + 864 MHC-I templates
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
