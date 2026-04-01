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
Module 8a: Immunogenicity + Foreignness Scoring
    BigMHC-IM (Nature Machine Intelligence 2023) → immunogenicity [0,1]
    k-mer foreignness vs UniProt human proteome (20,659 proteins)
    Agretopicity: log2(wt_rank / mut_rank)
    │
    ▼
Module 8b: Structural Scoring (Tiered)
    Tier 1: TCR-facing position lookup (P3-P7 for 9-mers) — instant
    Tier 2: PANDORA homology modeling + SASA — minutes, CPU
    Tier 3: AlphaFold via Tamarind.bio — hours, top 5-10 only
    │
    ▼
Module 9: Composite Ranking & Selection
    Weighted formula: immunogenicity (0.35) + foreignness (0.15)
      + agretopicity (0.10) + binding (0.15) + stability (0.10)
      + expression (0.10) + CCF (0.05)
    Bonuses: frameshift (+0.10), shared neoantigen (+0.15), CD4 (+0.05)
    Hard filters → top 20 with HLA diversity
    │
    ▼
Module 10: Polyepitope Design
    Greedy TSP ordering (minimize junctional epitopes)
    AAY linkers (MHC-I) / GPGPG linkers (MHC-II)
    Signal peptide + MITD domain
    │
    ▼
Module 11: mRNA Design
    LinearDesign: joint codon + mRNA structure optimization
    5' UTR (HBB) + CDS + 3' UTR (AES/mtRNR1) + poly-A (120nt)
    N1-methylpseudouridine substitution + CleanCap-AG
    │
    ▼
Outputs: mrna_sequence.fasta + synthesis_spec.json
```

## Tools

| Step | Tool | What it does | License |
|------|------|-------------|---------|
| Protein sequences | pyensembl (Ensembl GRCh38) | Real protein sequences from gene + mutation annotation | Apache 2.0 |
| MHC-I binding | MHCflurry 2.2 | Binding affinity + antigen processing + presentation | Apache 2.0 |
| Immunogenicity | BigMHC-IM | Transformer-based immunogenicity prediction | Academic |
| Foreignness | k-mer vs UniProt proteome | Self-similarity scoring (10.4M 9-mers) | -- |
| Structural | Position lookup / PANDORA / AlphaFold | TCR exposure scoring (tiered) | Mixed |
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
4. Runs MHCflurry on wildtype peptides at the same positions for agretopicity calculation

### Module 8a: Immunogenicity Scoring

Three independent signals for each peptide-HLA pair:

- **BigMHC-IM** (Albert et al., Nature Machine Intelligence 2023): Transformer-based immunogenicity prediction. Ensemble of 7 models trained on eluted ligand + immunogenicity data. Supports 8-14mer peptides with fuzzy HLA allele matching.

- **Foreignness**: k-mer overlap against the full UniProt human reference proteome (20,659 proteins, 10.4M 9-mers). Checks exact and 1-mismatch matches. Score of 0 = exact self-match, 1 = maximally foreign.

- **Agretopicity**: `log2(wt_binding_rank / mut_binding_rank)`. High positive = mutation creates a new binding event. Negative = wildtype already binds (tolerance risk).

### Module 8b: Structural Scoring

Tiered approach to score whether the mutated residue is TCR-exposed:

- **Tier 1** (instant): Position lookup table based on crystal structure surveys. TCR-facing positions are P4-P7 for 9-mers, with weighted scores reflecting average solvent exposure.

- **Tier 2** (1-2 min, CPU): PANDORA homology modeling generates a pMHC structure, then BioPython ShrakeRupley computes per-residue solvent-accessible surface area at the mutation site. Requires MODELLER license (free academic).

- **Tier 3** (hours, paid): AlphaFold batch prediction via Tamarind.bio API for top 5-10 candidates only. Full pMHC complex structure with real SASA scoring.

### Module 9: Ranking & Selection

Composite scoring formula with configurable weight profiles (`high_tmb`, `low_tmb`, `research`). Hard filters remove candidates with weak binding, low expression, or self-similarity. Top 20 selected with HLA allele diversity enforcement.

### Module 10: Polyepitope Design

Assembles selected neoantigens into an optimized polyepitope construct:
- Greedy TSP ordering to minimize junctional neoepitope formation
- AAY linkers (MHC-I) / GPGPG linkers (MHC-II)
- Signal peptide for secretory pathway + MITD domain for MHC-I loading
- Junctional QC flags high-risk junctions

### Module 11: mRNA Design

Converts the polyepitope protein into a synthesis-ready mRNA using **LinearDesign** for joint optimization of codon usage (CAI) and mRNA secondary structure (MFE). Adds 5' UTR, 3' UTR, poly-A tail, and outputs a complete synthesis specification with pseudouridine positions.

## Example Output (TCGA-EE-A3J5 Melanoma)

```
Patient: TCGA-EE-A3J5 (TCGA SKCM melanoma)
Variants: 2,584 somatic → 1,890 peptide windows → 260 strong binders
  → 84 pass filters → 20 selected
Polyepitope: 322 amino acids
mRNA: 1,253 nt | GC: 57.2% | MFE: -679.9 kcal/mol | CAI: 0.900
```

## Project Structure

```
ChatNAV/
├── bin/                            # Pipeline scripts
│   ├── maf_to_pipeline_input.py        # Step 0: MAF → peptides + MHCflurry
│   ├── score_immunogenicity.py         # Module 8a: BigMHC-IM + foreignness
│   ├── structural_scoring.py           # Module 8b: Tiered structural scoring
│   ├── rank_and_select.py              # Module 9: Composite ranking
│   ├── design_polyepitope.py           # Module 10: Polyepitope assembly
│   └── design_mrna.py                  # Module 11: mRNA + LinearDesign
├── conf/                           # Configuration
│   └── scoring_weights.yaml            # Weight profiles
├── reference/                      # Small reference data
│   ├── signal_peptides.fasta
│   ├── shared_neoantigens/
│   └── utr_library/
├── LinearDesign/                   # Submodule: codon optimization
├── models/DeepImmuno/              # Submodule: legacy immunogenicity model
├── docker/                         # Dockerfiles per module
├── modules/                        # Nextflow DSL2 modules
├── subworkflows/                   # Nextflow subworkflows
├── tests/
│   ├── run_e2e.py                      # Full e2e test
│   ├── test_scoring.py
│   ├── test_ranking.py
│   ├── test_polyepitope.py
│   └── test_mrna.py
├── main.nf                         # Nextflow entry point
├── nextflow.config
├── SETUP.md                        # Detailed installation log
└── plan.md                         # System design document
```

## Benchmark Results

Benchmarked against Gartner et al. 2021 (Nature Cancer) ground truth — 120 experimentally confirmed immunogenic peptide-HLA pairs vs 119 matched wildtype negatives from 61 NCI patients.

| Scoring Method | AUC | PPV@5 | PPV@10 | PPV@20 | PPV@50 |
|---|---|---|---|---|---|
| BigMHC-IM alone | 0.619 | 0.600 | 0.600 | 0.550 | 0.660 |
| Structural (Tier 1 position) | 0.883 | 1.000 | 1.000 | 1.000 | 0.980 |
| Composite (0.7 IM + 0.3 Struct) | 0.742 | 1.000 | 1.000 | 0.850 | 0.740 |

The structural tier 1 score is strong on this benchmark because confirmed immunogenic mutations are disproportionately at TCR-facing positions (P4-P7), while wildtype peptides score 0 by definition. A harder benchmark against all screened mutant candidates (including non-immunogenic mutants) is in progress.

## Known Limitations

- **CCF = 1.0**: Cancer cell fraction requires PyClone-VI with BAM files (VAF + CNV data). From MAF-only input, all mutations are assumed clonal.
- **No MHC Class II typing**: OptiType supports Class I only. Class II HLA typing requires HLA-HD (academic license).
- **MHCflurry vs NetMHCpan**: Uses MHCflurry (Apache 2.0). The 2% binding rank threshold is calibrated against MHCflurry's own percentile distribution.
- **LinearDesign on macOS ARM**: Requires Docker (linux/amd64) due to pre-compiled x86_64 shared libraries.

## References

- Albert et al. (2023). BigMHC: accurate prediction of MHC-I antigen presentation and immunogenicity. *Nature Machine Intelligence*.
- O'Donnell et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides. *Cell Systems*.
- Zhang et al. (2023). LinearDesign: Efficient algorithms for optimized mRNA sequence design. *Nature*.
- Gartner et al. (2021). A machine learning model for ranking candidate neoantigens. *Nature Medicine*.

## License

Pipeline code: MIT. External tools (BigMHC, LinearDesign, MHCflurry) have their own licenses -- see SETUP.md.
