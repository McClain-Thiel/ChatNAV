# Neoantigen Vaccine Pipeline — System Design Document

**Version:** 0.1.0  
**Status:** Draft  
**Date:** 2026-03-31

---

## Table of Contents

1. [Overview](#1-overview)
2. [Background & Scientific Rationale](#2-background--scientific-rationale)
3. [Pipeline Architecture](#3-pipeline-architecture)
4. [Module Specifications](#4-module-specifications)
5. [Infrastructure](#5-infrastructure)
6. [Data Design](#6-data-design)
7. [Caching Strategy](#7-caching-strategy)
8. [API Design](#8-api-design)
9. [Demo Data](#9-demo-data)
10. [Cost Model](#10-cost-model)
11. [Future Work](#11-future-work)

---

## 1. Overview

This document describes the architecture for a modular, cloud-native pipeline that takes genomic sequencing data from a cancer patient and produces a personalized mRNA neoantigen vaccine sequence ready for synthesis.

The pipeline accepts inputs ranging from raw FASTQ reads down to pre-called VCFs, enabling flexible entry points for different clinical and research contexts. Each module is independently swappable — as the science improves (particularly immunogenicity prediction), components can be upgraded without rebuilding the system.

The end-to-end output is a nucleotide string encoding a polyepitope mRNA construct, along with a synthesis specification sheet ready to hand to a manufacturing facility.

### Key design goals

- **Modular** — clean interfaces between steps, hot-swappable components
- **Resumable** — module-level caching, never re-run expensive steps unnecessarily
- **Flexible entry points** — accept FASTQ, BAM, VCF, or pre-processed inputs
- **Demo-ready** — pre-cached TCGA demo patients that run in <15 minutes
- **Cost-efficient** — ~$60-70 per patient end-to-end on AWS

---

## 2. Background & Scientific Rationale

### 2.1 What is a neoantigen vaccine?

Tumor cells accumulate somatic mutations not present in normal tissue. Some of these mutations produce altered proteins. If the immune system can be trained to recognize these altered proteins as foreign, it will selectively attack tumor cells displaying them.

A personalized neoantigen vaccine:
1. Identifies tumor-specific mutations via sequencing
2. Predicts which mutant peptides (neoantigens) will be displayed on the cell surface via MHC molecules
3. Encodes those peptides in an mRNA construct
4. Delivers the mRNA via lipid nanoparticles (LNPs) — the same technology as COVID vaccines
5. Trains the patient's immune system to hunt cells displaying those peptides

Unlike targeted therapies or checkpoint inhibitors, this approach is patient-specific by design: every vaccine is unique to the individual's tumor mutational landscape and HLA alleles.

### 2.2 Clinical evidence

The clearest clinical validation comes from the Moderna/Merck Phase IIb KEYNOTE-942 trial (mRNA-4157/intismeran autogene), which showed a **49% reduction in melanoma recurrence or death** over five years when a personalized neoantigen vaccine was combined with pembrolizumab (Keytruda). This used the same mRNA+LNP platform described here.

### 2.3 Why the pipeline is buildable now

Three converging technologies make this pipeline tractable today:

| Technology | Status |
|---|---|
| Whole genome sequencing | ~$500-1000/sample, 1-2 week turnaround |
| MHC binding prediction (MHCflurry 2.0) | Mature, well-validated, fully open source (Apache 2.0) |
| mRNA synthesis + LNP delivery | Industrialized post-COVID |

The weak link is **immunogenicity prediction** — predicting which MHC-binding peptides will actually activate T cells. Current SOTA (DeepImmuno, PRIME2.0, BigMHC) achieves ~PPV 0.37 on curated benchmarks. This is workable in practice because we encode 15-20 neoantigens per vaccine: even at PPV 0.37, the probability of zero immunogenic hits in 20 candidates is <0.002%.

The immunogenicity scoring module (Module 8) is explicitly designed as the hot-swappable research surface of this system. The current implementation uses BigMHC-IM (ensemble of 7 transformer models) as the primary immunogenicity predictor, replacing the earlier DeepImmuno-CNN. A LightGBM reranker (EXP-013) combines 13 features for candidate ranking.

### 2.4 Key scientific decisions

**Clonality filtering:** We target truncal mutations (present in all tumor cells, CCF > 0.8) rather than subclonal mutations. Targeting subclonal mutations allows the bulk of the tumor to escape vaccination. CCF is estimated from copy-number-corrected VAF using PyClone-VI on single biopsies. Multi-region sampling remains the gold standard but is operationally impractical for most clinical settings.

**Expanded antigen universe:** The pipeline searches beyond SNVs to include frameshift/indel neoantigens (no wildtype counterpart → no central tolerance), gene fusions (detectable from RNA-seq), and a shared neoantigen lookup table for recurrent driver mutations (KRAS G12V/D, TP53 R175H, IDH1 R132H) stratified by HLA allele.

**Structural scoring (optional):** AlphaFold2-Multimer predictions of pMHC complex structure can identify TCR-facing residue positions. Mutations at TCR-facing positions are more immunogenic than those buried in the MHC groove. This step is optional (adds cost and latency) and served via the NVIDIA API rather than local GPU.

---

## 3. Pipeline Architecture

### 3.1 Overview

```
INPUTS
  Tumor FASTQ / BAM / VCF
  Normal FASTQ / BAM
  RNA-seq FASTQ [optional]
  HLA alleles [optional, inferred if absent]
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 1+2: PREPROCESSING + VARIANT CALLING            │
│  Path A (CPU): BWA-MEM2 → GATK Mutect2    ~6-8hrs       │
│  Path B (GPU): Parabricks fq2bam+mutect   ~45min        │
│  Starting from BAM: skip to module 2                    │
│  Starting from VCF: skip to module 3                    │
└─────────────────────────┬───────────────────────────────┘
                          │  Annotated VCF, CNV, fusions
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 3: HLA TYPING                                   │
│  OptiType on BAM/FASTQ (Class I only)      ~30min       │
│  Skip if HLA alleles provided                           │
└──────────┬──────────────────────────────────────────────┘
           │  HLA alleles (A, B, C + Class II)
           │
           │   ┌──────────────────────────────────────────┐
           │   │  MODULE 4: EXPRESSION QUANTIFICATION     │
           │   │  Salmon on RNA-seq FASTQ       ~20min    │
           │   │  Fallback: GTEx priors by tumor type     │
           │   └──────────────────┬───────────────────────┘
           │                      │  TPM matrix
           ▼                      ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 5: CLONALITY ESTIMATION                         │
│  PyClone-VI → CCF per mutation             ~20min       │
│  Multi-sample: VCF intersection                         │
│  Shared neoantigen lookup (KRAS/TP53/IDH1)              │
└─────────────────────────┬───────────────────────────────┘
                          │  CCF-annotated VCF
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 6: NEOANTIGEN CANDIDATE GENERATION              │
│  pVACtools: SNV/indel peptides             ~20min       │
│  Frameshift novel ORFs                                  │
│  Fusion junction peptides [if RNA-seq]                  │
│  ERV peptides [optional]                                │
└─────────────────────────┬───────────────────────────────┘
                          │  FASTA of candidate 8-14mers
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 7: MHC BINDING PREDICTION                       │
│  MHCflurry 2.0 (MHC-I + processing)      ~10min        │
│  MHCnuggets (MHC-II)                                   │
│  All open-source, pip-installable                       │
└─────────────────────────┬───────────────────────────────┘
                          │  Binding ranks per (peptide, HLA)
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 8: IMMUNOGENICITY SCORING          ◄── RESEARCH │
│  DeepImmuno-CNN (open source)              ~15min        │
│  Foreignness score (IEDB similarity)                    │
│  Agretopicity (mutant vs wildtype binding delta)        │
│  Structural scoring via NVIDIA API [optional]           │
│    ESMFold pre-filter → AlphaFold2-Multimer top N       │
│    → TCR-facing residue exposure score                  │
└─────────────────────────┬───────────────────────────────┘
                          │  Immunogenicity scores per candidate
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 9: CANDIDATE RANKING + SELECTION                │
│  Weighted scoring across all signals      ~5min         │
│  Filter: MHC rank <2%, TPM >1, CCF >0.5               │
│  Prefer: CCF >0.8, frameshift, shared neoantigen        │
│  Output: Top 15-20 neoantigens                          │
└─────────────────────────┬───────────────────────────────┘
                          │  Ranked neoantigen list
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 10: POLYEPITOPE DESIGN                          │
│  Epitope ordering optimization (OptiVax)  ~5min         │
│  Linker selection (AAY/GPGPG)                           │
│  Junctional epitope check                               │
│  Signal peptide + MITD addition                         │
└─────────────────────────┬───────────────────────────────┘
                          │  Polyepitope amino acid sequence
                          ▼
┌─────────────────────────────────────────────────────────┐
│  MODULE 11: mRNA DESIGN                                 │
│  Codon optimization (LinearDesign)        ~5min         │
│  UTR selection (5' and 3')                              │
│  Pseudouridine substitution flag                        │
│  Poly-A tail specification                              │
└─────────────────────────┬───────────────────────────────┘
                          │
                          ▼
OUTPUTS
  mrna_sequence.fasta           Final nucleotide string
  final_neoantigens.tsv         Ranked table with all scores
  polyepitope.faa               Amino acid construct
  synthesis_spec.json           Spec sheet for synthesis facility
  qc_report.html                Per-module QC
  run_manifest.json             Full provenance
```

### 3.2 Entry points

The pipeline accepts inputs at four levels of pre-processing:

| Entry point | Modules skipped | Wall time | Use case |
|---|---|---|---|
| Raw FASTQ | None | ~2-3hrs (GPU) / ~7hrs (CPU) | Full clinical run |
| Aligned BAM | Module 1 | ~1.5-2hrs | Sequencing facility delivers BAM |
| Somatic VCF + expression | Modules 1-4 | ~45min | Pre-processed research data |
| Candidates FASTA + HLA | Modules 1-6 | ~20min | Demo / benchmarking |

Entry point is inferred automatically from which files are present in `inputs/`.

### 3.3 Module interface contracts

Each module consumes and produces files with fixed schemas. Swapping a module requires only that the replacement honors the same input/output schema.

```
Module 2 output:
  somatic_annotated.vcf.gz      standard VCF with VEP INFO fields
  cnv_segments.tsv              chr, start, end, log2_ratio, ploidy
  fusions.tsv                   gene1, gene2, junction_seq, tpm [if RNA-seq]

Module 3 output:
  hla_alleles.txt               one allele per line, e.g. HLA-A*02:01

Module 4 output:
  expression_tpm.tsv            gene_id, transcript_id, tpm, num_reads

Module 5 output:
  clonality_ccf.tsv             variant_id, ccf, cluster_id, is_truncal
                                + shared_neoantigen flag

Module 6 output:
  candidates.fasta              peptide sequences (8-14mers)
  candidates_meta.tsv           peptide_id, source_mutation, mutation_type,
                                gene, ccf, tpm, wildtype_peptide, is_frameshift,
                                is_shared, fusion_partner

Module 7 output:
  binding_predictions.tsv       peptide_id, hla_allele, mhc_class,
                                binding_rank, stability_rank,
                                cleavage_score, tap_score

Module 8 output:
  immunogenicity_scores.tsv     peptide_id, hla_allele,
                                immunogenicity_score, foreignness_score,
                                agretopicity, structural_score [optional]

Module 9 output:
  selected_neoantigens.tsv      top 15-20 neoantigens with all scores,
                                ordered by composite rank

Module 10 output:
  polyepitope.faa               FASTA of full polyepitope AA sequence
  polyepitope_design.tsv        epitope order, linker choices,
                                junctional epitope QC

Module 11 output:
  mrna_sequence.fasta           final nucleotide string
  synthesis_spec.json           codon table, modifications, UTR sequences,
                                poly-A length, pseudouridine positions
```

---

## 4. Module Specifications

### Module 1+2: Preprocessing + Variant Calling

**Purpose:** Convert raw FASTQ reads to an annotated somatic VCF.

**Path A — CPU (default)**
- Alignment: BWA-MEM2 against GRCh38
- Duplicate marking: GATK MarkDuplicates
- Base quality recalibration: GATK BQSR
- Somatic variant calling: GATK Mutect2 (scatter by chromosome, gather)
- CNV: GATK CNV pipeline or CNVkit
- Fusion detection: STAR-Fusion (if RNA-seq provided)
- Annotation: Ensembl VEP
- Instance: `r5.8xlarge` (32 vCPU, 256GB RAM)
- Wall time: ~6-8 hours
- Cost: ~$20

**Path B — GPU (Parabricks)**
- `parabricks fq2bam`: GPU-accelerated BWA + sorting + dedup + BQSR
- `parabricks mutectcaller`: GPU-accelerated Mutect2
- Requires NVIDIA GPU instance; free for research (<$5M revenue), ~$10K+/yr commercial license
- Instance: `p3.2xlarge` (1× V100)
- Wall time: ~45 minutes
- Cost: ~$8
- Selection: set `parabricks: true` in run config

**Output filters applied:**
- PASS variants only
- Allele frequency > 0.05
- Coverage > 10× at variant site
- Nonsynonymous consequence (VEP: missense, frameshift, stop_gained, splice_site, fusion)

---

### Module 3: HLA Typing

**Purpose:** Determine the patient's HLA alleles at 4-digit resolution.

- Tool: OptiType 1.3.5 (BSD-3 license, open source, `pip install` / Docker)
- Input: WGS/WES/RNA-seq BAM or FASTQ
- Alleles called: HLA-A, B, C (Class I only)
- **Limitation:** OptiType does not call Class II alleles (DRB1, DQB1). Class II HLA typing is not currently available in the pipeline. This means MHC-II binding predictions (Module 7, MHCnuggets) require Class II alleles to be provided externally via `hla_alleles.txt`. If Class II alleles are not provided, the pipeline runs with Class I only — the vaccine construct will lack CD4+ epitopes, which is acceptable but suboptimal. Adding arcasHLA (MIT license, supports partial Class II) or upgrading to HLA-HD (free academic license from Kyoto University) are future options.
- TCGA patients: use pre-published OptiType calls from GDC Pan-Immune study (bundled as reference)

---

### Module 4: Expression Quantification

**Purpose:** Estimate gene and transcript expression in the tumor.

- Tool: Salmon v1.10 (quasi-mapping, fast)
- Reference: GENCODE v45 transcriptome
- Output: TPM per gene and transcript
- Allele-specific expression: optional (requires phased VCF, adds complexity)
- **Fallback (no RNA-seq):** Use TCGA median TPM by tumor type from GTEx/TCGA combined reference. Flags all expression values as `imputed=true` in output. Candidates from imputed expression are down-weighted in Module 9 scoring.

---

### Module 5: Clonality Estimation

**Purpose:** Estimate the cancer cell fraction (CCF) of each mutation to prioritize truncal (clonal) neoantigens over subclonal ones.

**Single-sample (default):**
- Tool: PyClone-VI
- Input: CNV-corrected VAF per variant
- Output: CCF per variant, cluster assignments
- Truncal threshold: CCF > 0.8 (configurable)
- Minimum threshold: CCF > 0.5 (below this, variants are excluded)

**Multi-sample (optional):**
- If multiple tumor BAMs/VCFs provided (different regions or timepoints)
- Simple intersection: variants present in all samples = truncal
- More robust than single-sample VAF modeling
- Enabled by providing `additional_tumor_vcfs` in run config

**Shared neoantigen lookup:**
- After CCF estimation, annotate variants against the shared neoantigen reference table
- Table schema: `(gene, mutation, hla_allele, validated_immunogenic, source_study)`
- Current entries: KRAS G12V, G12D, G12R; TP53 R175H, R248W; IDH1 R132H; BRAF V600E
- Shared neoantigen flag receives a scoring bonus in Module 9 (configurable weight)

---

### Module 6: Neoantigen Candidate Generation

**Purpose:** Generate all candidate peptide sequences from the filtered, annotated, CCF-filtered variant list.

- Tool: pVACtools (pVACseq for SNV/indel, pVACfuse for fusions)
- Peptide window: 8-14mers (MHC-I), 13-25mers (MHC-II)
- Mutation types processed:
  - SNV / MNV: standard flanking window ±10aa
  - In-frame indel: flanking window around insertion/deletion
  - Frameshift: novel ORF from frameshift position to next stop codon
  - Fusion junction: junction-spanning peptides from STAR-Fusion output
  - ERV: optional, requires ERV reference (hg38-ERV annotations)
- Wildtype peptide generated for each candidate (for agretopicity calculation)
- Peptides matching wildtype proteome exactly are discarded

**Why frameshifts matter:** Frameshift mutations produce completely novel protein sequences downstream of the frameshift — no wildtype counterpart exists. This means no central tolerance has deleted T cells that recognize them. Frameshift neoantigens are enriched for immunogenicity per mutation relative to SNVs and receive a scoring bonus in Module 9.

---

### Module 7: MHC Binding Prediction

**Purpose:** Predict which candidate peptides will be presented on the cell surface via MHC molecules.

**All tools are open-source and pip-installable — no academic licenses required.**

**MHC-I (CD8+ T cell pathway):**
- Primary: MHCflurry 2.0 (Apache 2.0 license, `pip install mhcflurry`)
- Pan-allele model supporting 14,000+ MHC-I alleles
- Integrates binding affinity + antigen processing (cleavage + TAP transport) in a single model
- Replaces four separate NetMHC tools: NetMHCpan, NetMHCstabpan, NetChop, NetCTLpan
- Benchmarks: comparable or better than NetMHCpan 4.0 on mass-spec ligand data; 396× faster
- Prediction speed: >7,000 predictions/second

**MHC-II (CD4+ T cell pathway):**
- MHCnuggets (MIT license, `pip install mhcnuggets`)
- Supports both Class I and Class II — used here for Class II only
- Requires Class II HLA alleles (DRB1, DQB1) in `hla_alleles.txt`
- If no Class II alleles provided, MHC-II prediction is skipped (vaccine is Class I only)
- CD4+ response amplifies and sustains the CD8+ response — including both pathways improves overall vaccine efficacy

**Filtering thresholds (configurable):**
```
MHC-I presentation percentile:  < 2%  (weak binder)
                                 < 0.5% (strong binder, preferred)
MHC-II binding rank:             < 2%
```

**Scale:** A typical run produces 50-200 MHC binders from ~500-2000 candidate peptides. All are carried forward to Module 8 for scoring.

---

### Module 8: Immunogenicity Scoring

**Purpose:** Rank MHC-binding candidates by predicted likelihood of activating a T cell response.

This is the primary research surface of the pipeline. The scoring is a composite of multiple signals.

**Signal 1 — DeepImmuno-CNN (primary immunogenicity model)**
- Open source (GitHub: frankligy/DeepImmuno), no license restrictions
- CNN-based model predicting CD8+ T cell immunogenicity for 9-10mer peptides
- Trained on experimentally validated immunogenic/non-immunogenic epitopes
- Outperforms IEDB immunogenicity predictor and DeepHLApan on benchmarks
- Input: peptide sequence (one-hot encoded)
- Output: immunogenicity probability score [0, 1]
- Fallback: biophysical heuristic based on TCR-facing residue properties (if model weights unavailable)

**Signal 2 — Foreignness score**
- IEDB similarity search: how different is the mutant peptide from anything in the normal human proteome?
- High foreignness → no central tolerance → higher predicted immunogenicity
- Computed as: 1 - max(sequence_similarity to self-proteome)

**Signal 3 — Agretopicity**
- Delta between mutant and wildtype MHC binding rank
- High agretopicity: mutation creates a new binding event (not just modifying existing one)
- Low agretopicity: wildtype also binds well → T cells may be tolerized
- Filter: mutant rank must be meaningfully better than wildtype rank

**Signal 4 — Structural scoring via NVIDIA API (optional)**
- AlphaFold2-Multimer prediction of pMHC complex structure
- Endpoint: `https://health.api.nvidia.com/v1/biology/deepmind/alphafold2-multimer`
- Input: MHC heavy chain sequence + beta-2 microglobulin + peptide (as multimer)
- Analysis: identify TCR-facing residue positions; score mutation at TCR-facing vs groove-buried position
- Pre-filter with ESMFold (faster, cheaper) to select top 20-30 candidates before AlphaFold
- Cost: ~$0.05-0.10/prediction; ~$2-5 per patient run
- Wall time: ~1-2 hours (API calls, partially parallelizable)
- Enabled by setting `structural_scoring: true` in run config

```python
# NVIDIA API call sketch
payload = {
    "sequences": [
        {"sequence": hla_heavy_chain_aa},         # e.g. HLA-A*02:01 sequence
        {"sequence": B2M_AA},                      # beta-2 microglobulin (constant)
        {"sequence": peptide_9mer}                 # candidate neoantigen
    ]
}
response = requests.post(
    "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2-multimer",
    headers={"Authorization": f"Bearer {NVIDIA_API_KEY}"},
    json=payload
)
# returns PDB structure → parse TCR-facing residue exposure
```

**Swap note:** The entire Module 8 scoring logic is encapsulated behind a single interface. To upgrade to a better model (e.g., fine-tuned ESM-based immunogenicity predictor, or a model trained on proprietary trial data), replace the scoring function while keeping input/output schemas identical.

---

### Module 9: Candidate Ranking and Selection

**Purpose:** Integrate all signals into a final ranked list and select the top 15-20 neoantigens for the vaccine.

**Composite score formula (configurable weights):**

```
composite_score = (
    w1 × immunogenicity_score      # default: 0.35
  + w2 × foreignness_score         # default: 0.15
  + w3 × agretopicity              # default: 0.10
  + w4 × (1 - mhc_binding_rank)   # default: 0.15
  + w5 × stability_score           # default: 0.10
  + w6 × log10(tpm + 1) / 5       # default: 0.10  (normalized expression)
  + w7 × ccf                       # default: 0.05
  + w8 × structural_score          # default: 0.00 (0.10 if structural enabled)
)

# Bonuses (additive)
+ 0.10 if is_frameshift
+ 0.15 if is_shared_neoantigen (pre-validated)
+ 0.05 if mhc_class_ii_binder     (CD4+ epitope available)
```

**Hard filters (applied before scoring):**
```
mhc_binding_rank    < 0.02
tpm                 > 1.0
ccf                 > 0.5
agretopicity        > 0      (mutant binds better than wildtype)
not_self_match      = True
```

**Tumor type weight profiles (configurable presets):**

| Profile | Notes |
|---|---|
| `high_tmb` (melanoma, NSCLC) | Default weights, volume is available |
| `low_tmb` (PDAC, GBM, MSS-CRC) | Upweight frameshift bonus, shared neoantigen, structural |
| `research` | Include all candidates, no top-N cap, full score breakdown |

---

### Module 10: Polyepitope Design

**Purpose:** Assemble selected neoantigens into an optimized polyepitope construct.

- Epitope ordering: OptiVax (minimizes junctional epitope formation)
- Linker selection:
  - `AAY` between MHC-I epitopes (clean proteasomal cleavage)
  - `GPGPG` between MHC-II epitopes
  - `AAAAA` at construct termini
- Junctional epitope check: sliding window across all linker junctions, score with MHCflurry; flag any junction peptide with rank < 0.5% for manual review
- Signal peptide: standard secretory signal (METDTLLLWVLLLWVPGSTGD) prepended to route through secretory pathway for MHC-II loading
- MITD (MHC targeting domain): optional, appended to C-terminus to improve both MHC-I and MHC-II presentation

---

### Module 11: mRNA Design

**Purpose:** Convert the polyepitope amino acid sequence into a synthesis-ready mRNA nucleotide string.

- Codon optimization: LinearDesign (Baidu, 2023) — jointly optimizes codon usage and mRNA secondary structure for maximum stability and translation efficiency
- Pseudouridine substitution: all uridines flagged for N1-methylpseudouridine substitution at synthesis (reduces innate immune sensing of exogenous mRNA — the Karikó/Weissman Nobel insight)
- UTR selection: 5' and 3' UTRs from validated high-expression reference library
- Poly-A tail: 120nt poly-A (configurable)
- Output: nucleotide FASTA + synthesis specification JSON

**Synthesis spec JSON format:**
```json
{
  "mrna_sequence": "AUGCUG...",
  "length_nt": 2847,
  "codon_optimization": "LinearDesign",
  "modifications": {
    "pseudouridine": true,
    "cap_analog": "CleanCap-AG",
    "poly_a_length": 120
  },
  "utr_5prime": "GGGAAAUAAGAGAGAAAAGAAGAG...",
  "utr_3prime": "UGAUAAUAGGCUCGAGCAUCGAGC...",
  "neoantigen_positions": [
    {"rank": 1, "gene": "BRAF", "mutation": "V600E", 
     "nt_start": 102, "nt_end": 129}
  ],
  "synthesis_notes": "N1-methylpseudouridine substitution throughout"
}
```

---

## 5. Infrastructure

### 5.1 Architecture diagram

```
┌─────────────────────────────────────────────────────────┐
│                    USER / CLIENT                        │
│              Browser (React/Next.js)                    │
└────────────────────────┬────────────────────────────────┘
                         │ HTTPS
┌────────────────────────▼────────────────────────────────┐
│                   API LAYER                             │
│              FastAPI on EC2 t3.medium                   │
│              (or ECS Fargate for auto-scale)            │
└──────────┬───────────────┬──────────────────────────────┘
           │               │
  ┌────────▼──────┐  ┌─────▼──────────┐
  │   PostgreSQL  │  │   S3 Bucket    │
  │   RDS t3.med  │  │   (see §6.1)   │
  │   Job state   │  │   All files    │
  │   Results     │  └────────────────┘
  └───────────────┘
           │
┌──────────▼──────────────────────────────────────────────┐
│                  AWS BATCH                              │
│           Nextflow as orchestrator                      │
│                                                         │
│  CPU queue:  r5.8xlarge, c5.4xlarge (spot)             │
│  GPU queue:  p3.2xlarge (Parabricks, spot)             │
│  GPU queue:  g4dn.xlarge (ESMFold pre-filter)          │
│                                                         │
│  Each module = Nextflow process = Docker container      │
└──────────┬──────────────────────────────────────────────┘
           │ API calls
┌──────────▼──────────────────────────────────────────────┐
│              NVIDIA API (external)                      │
│   AlphaFold2-Multimer, ESMFold                         │
│   health.api.nvidia.com                                 │
└─────────────────────────────────────────────────────────┘
```

### 5.2 Nextflow pipeline structure

```groovy
// main.nf sketch

nextflow.enable.dsl=2

workflow {
    // Infer entry point from available inputs
    entry_point = detect_entry_point(params)

    if (entry_point <= 2) {
        // Modules 1+2: preprocessing + variant calling
        if (params.parabricks) {
            (bam_tumor, bam_normal) = parabricks_fq2bam(
                params.tumor_fastq, params.normal_fastq)
            vcf = parabricks_mutect(bam_tumor, bam_normal)
        } else {
            bam_tumor = align(params.tumor_fastq, params.reference)
            bam_normal = align(params.normal_fastq, params.reference)
            vcf = mutect2_scatter_gather(bam_tumor, bam_normal)
        }
        vcf_annotated = vep_annotate(vcf)
    } else {
        vcf_annotated = Channel.fromPath(params.vcf)
    }

    // Module 3: HLA typing (parallel with variant calling)
    if (entry_point <= 3 && !params.hla_alleles) {
        hla = hla_hd(bam_normal ?: params.normal_bam)
    } else {
        hla = Channel.fromPath(params.hla_alleles)
    }

    // Module 4: expression (parallel)
    if (params.rnaseq_fastq && entry_point <= 4) {
        expression = salmon_quant(params.rnaseq_fastq)
    } else {
        expression = gtex_fallback(params.tumor_type)
    }

    // Module 5: clonality
    ccf_vcf = pyclone_vi(vcf_annotated, expression)
           | shared_neoantigen_lookup

    // Module 6: candidate generation
    candidates = pvactools(ccf_vcf, hla)

    // Module 7: MHC binding
    binding = netmhcpan(candidates, hla)
            | netmhcstabpan(hla)
            | netchop
            | netctlpan

    // Module 8: immunogenicity
    scores = prime(binding)
           | foreignness_score
           | agretopicity(binding)

    if (params.structural_scoring) {
        scores = esmfold_prefilter(scores)
               | alphafold_multimer_api   // NVIDIA API
               | tcr_exposure_score
    }

    // Modules 9-11: ranking + design
    selected = rank_and_filter(scores)
    polyepitope = optiVax_design(selected)
    mrna = linear_design(polyepitope)
         | generate_synthesis_spec

    publish_outputs(mrna, selected, polyepitope)
}
```

### 5.3 Instance types and AWS Batch queues

| Queue | Instance type | Use | Spot? |
|---|---|---|---|
| `cpu-standard` | `r5.8xlarge` | BWA-MEM2 alignment | Yes |
| `cpu-scatter` | `c5.4xlarge` × 24 | Mutect2 scatter | Yes |
| `gpu-parabricks` | `p3.2xlarge` | Parabricks | Yes |
| `gpu-inference` | `g4dn.xlarge` | ESMFold pre-filter | Yes |
| `cpu-light` | `c5.2xlarge` | All other modules | Yes |

All queues use Spot instances with on-demand fallback. Spot interruption handling is managed by Nextflow's `-resume` capability — interrupted jobs restart from the last completed process.

### 5.4 Docker containers (one per module)

```
neoantigen/preprocessing:0.1.0      BWA-MEM2, samtools, GATK
neoantigen/parabricks:0.1.0         Parabricks 4.1 (requires NVIDIA base)
neoantigen/optitype:0.1.0            OptiType 1.3.5 (BSD-3, Class I only)
neoantigen/salmon:0.1.0             Salmon 1.10
neoantigen/pyclone:0.1.0            PyClone-VI
neoantigen/pvactools:0.1.0          pVACtools 4.x
neoantigen/mhc-binding:0.1.0         MHCflurry 2.0 (Apache-2.0) + MHCnuggets (MIT)
neoantigen/scoring:0.1.0            DeepImmuno-CNN, foreignness, agretopicity
neoantigen/structural:0.1.0         ESMFold + NVIDIA API client
neoantigen/design:0.1.0             OptiVax, LinearDesign, synthesis spec
```

Each container is versioned and immutable — upgrading a tool means bumping the container version and invalidating downstream caches.

---

## 6. Data Design

### 6.1 S3 layout

```
s3://{bucket}/

  reference/                          # static reference data, ~500GB total
    hg38/
      GRCh38.fa.gz + indices
    vep_cache/                        # Ensembl VEP cache v110, ~50GB
    netmhcpan/                        # allele models
    gencode_v45/                      # transcriptome for Salmon
    shared_neoantigens/
      kras_hla_lookup.tsv
      tp53_hla_lookup.tsv
      idh1_hla_lookup.tsv
    tcga/
      hla_types_shukla2015.tsv        # pre-published TCGA HLA types
      skcm_expression/                # pre-downloaded TCGA-SKCM TPM
      skcm_vcfs/                      # pre-downloaded masked somatic VCFs

  patients/
    {patient_id}/
      meta.json
      inputs/                         # user-supplied files
        tumor.fastq.gz
        normal.fastq.gz
        rnaseq.fastq.gz
        tumor.bam
        normal.bam
        somatic.vcf.gz
        expression.tsv
        hla_alleles.txt
      cache/
        module_02/
          somatic_annotated.vcf.gz + .tbi
          cnv_segments.tsv
          fusions.tsv
        module_03/
          hla_alleles.txt
        module_04/
          expression_tpm.tsv
        module_05/
          clonality_ccf.tsv
        module_06/
          candidates.fasta
          candidates_meta.tsv
        module_07/
          binding_predictions.tsv
        module_08/
          immunogenicity_scores.tsv
          structural_scores.tsv
      outputs/
        final_neoantigens.tsv
        polyepitope.faa
        mrna_sequence.fasta
        synthesis_spec.json
        qc_report.html
        run_manifest.json

  demo/
    tcga_skcm/
      TCGA-DF-A2KN/                   # melanoma, high TMB
      TCGA-EE-A3J5/
      TCGA-ER-A19O/
      TCGA-GN-A262/
    gartner_nci/                      # with ground truth immunogenicity labels
      NCI_patient_001/
      NCI_patient_002/
      NCI_patient_003/
```

### 6.2 PostgreSQL schema

```sql
-- Patients
CREATE TABLE patients (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    external_id     TEXT UNIQUE,
    tumor_type      TEXT,
    created_at      TIMESTAMPTZ DEFAULT NOW(),
    metadata        JSONB
);

-- Pipeline runs
CREATE TABLE pipeline_runs (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    patient_id      UUID REFERENCES patients(id),
    created_at      TIMESTAMPTZ DEFAULT NOW(),
    completed_at    TIMESTAMPTZ,
    status          TEXT CHECK (status IN (
                        'queued','running','complete','failed')),
    config          JSONB,
    starting_module INT DEFAULT 1,
    parabricks      BOOLEAN DEFAULT FALSE,
    structural      BOOLEAN DEFAULT FALSE,
    error           TEXT,
    s3_prefix       TEXT
);

-- Per-module tracking
CREATE TABLE module_runs (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    run_id          UUID REFERENCES pipeline_runs(id),
    module_number   INT,
    module_name     TEXT,
    status          TEXT,
    started_at      TIMESTAMPTZ,
    completed_at    TIMESTAMPTZ,
    batch_job_id    TEXT,
    instance_type   TEXT,
    cost_usd        NUMERIC(10,4),
    cache_hit       BOOLEAN DEFAULT FALSE,
    output_s3_key   TEXT,
    error           TEXT
);

-- Neoantigen results
CREATE TABLE neoantigens (
    id                      UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    run_id                  UUID REFERENCES pipeline_runs(id),
    patient_id              UUID REFERENCES patients(id),
    rank                    INT,
    peptide_sequence        TEXT NOT NULL,
    source_gene             TEXT,
    mutation                TEXT,
    mutation_type           TEXT CHECK (mutation_type IN (
                                'snv','indel','frameshift','fusion','erv')),
    hla_allele              TEXT,
    mhc_class               TEXT CHECK (mhc_class IN ('I','II')),
    mhc_binding_rank        NUMERIC(5,3),
    mhc_stability_rank      NUMERIC(5,3),
    expression_tpm          NUMERIC(10,3),
    expression_imputed      BOOLEAN DEFAULT FALSE,
    ccf                     NUMERIC(5,3),
    immunogenicity_score    NUMERIC(8,5),
    foreignness_score       NUMERIC(8,5),
    agretopicity            NUMERIC(8,5),
    structural_score        NUMERIC(8,5),
    composite_score         NUMERIC(8,5),
    is_shared_neoantigen    BOOLEAN DEFAULT FALSE,
    is_frameshift           BOOLEAN DEFAULT FALSE,
    is_cd4_epitope          BOOLEAN DEFAULT FALSE,
    selected                BOOLEAN DEFAULT FALSE,
    wildtype_peptide        TEXT,
    scores_detail           JSONB
);

-- Final outputs
CREATE TABLE pipeline_outputs (
    id                  UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    run_id              UUID REFERENCES pipeline_runs(id),
    polyepitope_aa      TEXT,
    mrna_sequence       TEXT,
    mrna_length_nt      INT,
    neoantigen_count    INT,
    synthesis_spec      JSONB,
    created_at          TIMESTAMPTZ DEFAULT NOW()
);

-- Module cache
CREATE TABLE module_cache (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    patient_id      UUID REFERENCES patients(id),
    module_number   INT,
    cache_key       TEXT,
    s3_key          TEXT,
    tool_version    TEXT,
    created_at      TIMESTAMPTZ DEFAULT NOW(),
    is_valid        BOOLEAN DEFAULT TRUE,
    UNIQUE(patient_id, module_number, cache_key)
);

-- Demo patients
CREATE TABLE demo_patients (
    id                          UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    patient_id                  UUID REFERENCES patients(id),
    source                      TEXT CHECK (source IN (
                                    'tcga_skcm','gartner_nci')),
    tcga_barcode                TEXT,
    precomputed_through_module  INT,
    ground_truth_immunogenic    JSONB,
    description                 TEXT
);

CREATE INDEX idx_neoantigens_run ON neoantigens(run_id);
CREATE INDEX idx_neoantigens_patient ON neoantigens(patient_id);
CREATE INDEX idx_module_runs_run ON module_runs(run_id);
CREATE INDEX idx_module_cache_lookup ON module_cache(patient_id, module_number);
```

---

## 7. Caching Strategy

### 7.1 Cache key generation

Each module's cache key is a deterministic hash of its inputs:

```python
def compute_cache_key(
    input_s3_etags: list[str],   # S3 ETag of each input file
    tool_version: str,            # e.g. "netmhcpan-4.1"
    config_params: dict           # module-specific params
) -> str:
    content = "|".join(sorted(input_s3_etags))
    content += f"|{tool_version}"
    content += f"|{json.dumps(config_params, sort_keys=True)}"
    return hashlib.sha256(content.encode()).hexdigest()
```

### 7.2 Cache invalidation rules

| Trigger | Action |
|---|---|
| Tool version bump | Invalidate that module + all downstream modules |
| Config param change | Invalidate affected module + downstream |
| New input file (different ETag) | Invalidate all modules |
| Manual invalidation | Set `is_valid = false` in module_cache |

### 7.3 Demo patient cache behavior

Demo patients have Modules 1-7 pre-computed and cached. When a demo run is submitted:

```
Module 1:  ✓ Cached  (preprocessing)          0s
Module 2:  ✓ Cached  (variant calling)         0s
Module 3:  ✓ Cached  (HLA typing)              0s
Module 4:  ✓ Cached  (expression)              0s
Module 5:  ✓ Cached  (clonality)               0s
Module 6:  ✓ Cached  (candidate generation)    0s
Module 7:  ✓ Cached  (MHC binding)             0s
Module 8:  ⟳ Running (immunogenicity scoring)  ~5min
Module 9:  ⟳ Queued  (ranking)
Module 10: ⟳ Queued  (polyepitope design)
Module 11: ⟳ Queued  (mRNA design)

Total demo time: ~10-15 minutes
```

---

## 8. API Design

### 8.1 Endpoints

```
POST   /patients
       Create patient record, get presigned S3 upload URLs

POST   /patients/{id}/runs
       Submit pipeline run with config

GET    /patients/{id}/runs/{run_id}/status
       Poll job status, per-module progress

GET    /patients/{id}/runs/{run_id}/results
       Retrieve ranked neoantigens + outputs

GET    /patients/{id}/runs/{run_id}/results/download
       Signed S3 URLs for output files

GET    /demo/patients
       List available demo patients with descriptions

POST   /demo/patients/{demo_id}/runs
       Submit demo run (uses cached modules 1-7)
```

### 8.2 Request / response examples

**Submit a run:**
```json
POST /patients/{id}/runs
{
  "starting_module": 3,
  "config": {
    "tumor_type": "SKCM",
    "parabricks": false,
    "structural_scoring": false,
    "top_n_neoantigens": 20,
    "weight_profile": "high_tmb",
    "min_ccf": 0.5,
    "min_tpm": 1.0,
    "binding_rank_threshold": 0.02
  }
}

→ { "run_id": "uuid", "estimated_completion": "2026-03-31T15:30:00Z" }
```

**Poll status:**
```json
GET /patients/{id}/runs/{run_id}/status
→ {
    "status": "running",
    "current_module": 8,
    "modules": {
      "1": {"status": "complete", "cache_hit": true,  "cost_usd": 0},
      "2": {"status": "complete", "cache_hit": true,  "cost_usd": 0},
      "7": {"status": "complete", "cache_hit": false, "cost_usd": 2.14},
      "8": {"status": "running",  "started_at": "2026-03-31T15:10:00Z"}
    },
    "estimated_completion": "2026-03-31T15:20:00Z"
  }
```

**Results:**
```json
GET /patients/{id}/runs/{run_id}/results
→ {
    "run_id": "uuid",
    "status": "complete",
    "neoantigen_count": 20,
    "neoantigens": [
      {
        "rank": 1,
        "peptide": "YLQPRTFLL",
        "gene": "BRAF",
        "mutation": "V600E",
        "mutation_type": "snv",
        "hla_allele": "HLA-A*02:01",
        "mhc_binding_rank": 0.003,
        "ccf": 0.94,
        "tpm": 847.2,
        "immunogenicity_score": 0.821,
        "composite_score": 0.743,
        "is_shared_neoantigen": true
      }
    ],
    "mrna_length_nt": 2847,
    "total_cost_usd": 58.40
  }
```

---

## 9. Demo Data

### 9.1 TCGA-SKCM demo patients (Option B)

Four melanoma patients from TCGA selected for high TMB and clear neoantigen signal. All open-access data from GDC portal — no dbGaP approval required.

| TCGA Barcode | TMB | Known neoantigens | Notes |
|---|---|---|---|
| TCGA-DF-A2KN | High | Yes | Frequently cited in neoantigen literature |
| TCGA-EE-A3J5 | High | Yes | Strong BRAF V600E |
| TCGA-ER-A19O | High | Yes | Multiple frameshift neoantigens |
| TCGA-GN-A262 | High | Yes | Good MHC-II epitopes |

**Data sourced from:**
- Masked somatic VCFs: GDC Data Portal, Data Type = "Masked Somatic Mutation"
- Expression (TPM): GDC Data Portal, Data Type = "Gene Expression Quantification"  
- HLA types: Shukla et al. 2015 (Nature Methods), supplementary table

Pre-computed through Module 7. Modules 8-11 run live on demo submission.

### 9.2 Gartner NCI dataset (Option A — validation)

112 patients with experimentally validated neoantigen immunogenicity from Gartner et al. 2022 (Nature). Ground truth T cell response data (ELISpot/IFNγ) available for each mutation.

- Source: Supplementary tables, Gartner et al. Nature 2022
- Use case: benchmark and validate Module 8 immunogenicity scoring
- Three patients selected for demo (pre-processed, ground truth labels available)
- Enables head-to-head comparison: pipeline prediction vs. experimentally measured immunogenicity

### 9.3 Downloading demo data

```bash
# TCGA-SKCM via GDC API (no auth required for masked somatic + expression)
gdc-client download \
  --manifest tcga_skcm_manifest.txt \
  --dir ./demo/tcga_skcm/

# Gartner NCI — from paper supplementary tables
# Gartner et al. Nature 2022: https://doi.org/10.1038/s41586-022-04674-5
# Supplementary Table 1: patient mutations + HLA types
# Supplementary Table 2: immunogenicity assay results
```

---

## 10. Cost Model

### 10.1 Per-patient compute costs

| Entry point | Path | Compute cost | Wall time |
|---|---|---|---|
| Raw FASTQ | CPU (BWA + Mutect2) | ~$28 | ~7 hours |
| Raw FASTQ | GPU (Parabricks) | ~$12 | ~2.5 hours |
| BAM provided | CPU | ~$14 | ~3.5 hours |
| VCF + expression | Any | ~$5 | ~45 min |
| Demo (cached 1-7) | Any | ~$2 | ~12 min |

Full run with structural scoring (AlphaFold NVIDIA API):
```
Base compute:           ~$12   (Parabricks path)
NVIDIA API calls:       ~$4    (50 pMHC predictions)
Storage I/O:            ~$2
Overhead (API, DB):     ~$2
─────────────────────────────
Total:                  ~$20 per patient
```

### 10.2 Storage costs

| Data type | Size | Monthly cost |
|---|---|---|
| Raw FASTQs (3 samples) | ~300GB | ~$7 |
| BAMs | ~200GB | ~$5 |
| Intermediates | ~10GB | ~$0.25 |
| Outputs | <1GB | negligible |
| Reference data (shared) | ~500GB | ~$12 (amortized) |

Recommendation: delete BAMs after pipeline completes, retain FASTQs and outputs. Total per-patient ongoing storage: ~$0.50/month.

### 10.3 Infrastructure fixed costs

| Resource | Monthly cost |
|---|---|
| EC2 t3.medium (API server) | ~$30 |
| RDS t3.medium (Postgres) | ~$50 |
| S3 reference data | ~$12 |
| Data transfer | ~$5-20 |
| **Total fixed** | **~$100-110/month** |

---

## 11. Future Work

### 11.1 Module 8 improvements (research surface)

- **Structural scoring default-on:** As NVIDIA API costs decrease and AlphaFold-Multimer call times improve, enable structural scoring by default. TCR-facing residue exposure is the most principled addition to immunogenicity prediction available today.
- **Fine-tuned language model:** Train an ESM2-based immunogenicity predictor on the combined NCI + TESLA + HiTIDE public datasets. Small dataset but a fine-tuned LM should outperform DeepImmuno-CNN's architecture.
- **PRIME 2.1 integration:** PRIME (academic license, Gfeller lab) can be added as an additional scoring signal alongside DeepImmuno for ensemble predictions.
- **Proprietary data integration:** If clinical outcome data becomes available (measured T cell responses from vaccinated patients), retrain Module 8 on that signal. This is the moat that Gritstone/BioNTech/Moderna hold.

### 11.2 Clonality improvements

- **Multi-region sampling support:** Pipeline already supports multi-sample VCF intersection. Build UI to accept multiple tumor biopsies and automatically use TRACERx-style clonality inference.
- **Liquid biopsy / ctDNA:** Track clonal evolution longitudinally from blood draws. Enables vaccine update without re-biopsy.

### 11.3 Antigen universe expansion

- **ERV neoantigens:** Endogenous retrovirus reactivation creates highly foreign peptides. Requires ERV reference annotation layer. High immunogenicity per antigen but bioinformatics pipeline is less mature.
- **Non-canonical translation:** Upstream ORFs, ribosomal frameshifting, cryptic peptides. Requires ribosome profiling (Ribo-seq) data, currently experimental.
- **Post-translational modifications:** Phosphopeptides, glycopeptides as antigens. Very early stage.

### 11.4 Combination therapy support

- **Checkpoint inhibitor pairing:** Output optimal checkpoint inhibitor recommendation based on tumor type and PD-L1 expression (from RNA-seq Module 4).
- **TME scoring:** Add tumor microenvironment immunosuppression score (TIDE, ESTIMATE) to flag cases where T cell exclusion is likely — may indicate need for combination therapy beyond checkpoint inhibition.

### 11.5 Manufacturing interface

- Direct API integration with mRNA synthesis facilities (currently manual PDF/JSON handoff)
- Turnaround time tracking from synthesis order to delivery
- Automated QC comparison between ordered and received sequence

---

*End of document. Version 0.1.0.*
