# Infrastructure — ChatNAV

Current deployment architecture and the path toward a self-service platform.

## Current Architecture

ChatNAV is a **Nextflow DSL2 pipeline** that runs locally or on AWS Batch. Each pipeline stage is containerized and orchestrated by Nextflow, which handles dependency resolution, parallelism, and resume-on-failure.

```
                         ┌──────────────────────────────┐
                         │        Nextflow (main.nf)     │
                         │   orchestrates 11 modules     │
                         └──────────┬───────────────────┘
                                    │
         ┌──────────────────────────┼──────────────────────────┐
         ▼                          ▼                          ▼
┌─────────────────┐      ┌──────────────────┐      ┌──────────────────┐
│ CPU stages       │      │ GPU stages        │      │ External APIs    │
│ MHCflurry        │      │ BigMHC-IM         │      │ AlphaFold2 NIM   │
│ pyensembl        │      │ Parabricks (opt)  │      │ (structural T3)  │
│ LinearDesign     │      │                   │      │                  │
│ polyepitope/mRNA │      │                   │      │                  │
└─────────────────┘      └──────────────────┘      └──────────────────┘
```

### Entry points

The pipeline auto-detects the starting point from available inputs:

| Input provided | Entry point | Modules run |
|---|---|---|
| FASTQ (tumor + normal + RNA) | Raw sequencing | 1-11 (full pipeline) |
| BAM files | Post-alignment | 2-11 |
| VCF + expression | Post-calling | 5-11 |
| MAF + expression + HLA | Vaccine design only | 8-11 |
| Binding predictions + metadata | Scoring only | 8-11 |

Most development and benchmarking uses the MAF entry point (modules 8-11), which takes seconds to minutes. The full FASTQ pipeline takes 4-8 hours per patient on AWS Batch.

### Compute resources

| Resource | Instance | Purpose | Cost |
|---|---|---|---|
| **chatnav-gpu** | g6e.4xlarge (L40S 46GB) | BigMHC scoring, benchmarking | ~$1.86/hr |
| **g6-big2** | g6.4xlarge | General GPU work | ~$1.32/hr |
| **AWS Batch** (CPU) | m5.4xlarge pool | MHCflurry, alignment, variant calling | on-demand |
| **AWS Batch** (GPU) | g4dn.xlarge pool | Parabricks GPU variant calling | on-demand |
| **S3** | s3://chatnav-pipeline-data/ | ~620GB cached reference data, intermediates | ~$14/mo |

### Docker images

14 Dockerfiles in `docker/`, one per pipeline stage:

```
docker/
├── preprocessing/     # BWA-MEM2 + samtools + GATK
├── parabricks/        # NVIDIA Parabricks (GPU alignment + calling)
├── vep/               # Ensembl VEP annotation
├── hla-hd/            # HLA-HD typing
├── salmon/            # Salmon RNA quantification
├── pyclone/           # PyClone-VI clonality
├── pvactools/         # pVACseq candidate generation
├── star-fusion/       # STAR-Fusion
├── cnvkit/            # CNVkit copy number
├── mhc-binding/       # MHCflurry + MHCnuggets
├── netmhc/            # NetMHCpan (academic license)
├── scoring/           # BigMHC-IM + foreignness + agretopicity
├── pandora/           # PANDORA homology modeling (structural T2)
└── design/            # LinearDesign + polyepitope + mRNA design
```

### Nextflow profiles

| Profile | Config | Use case |
|---|---|---|
| `local` | `conf/local.config` | Development, single patient, macOS/Linux |
| `aws` | `conf/aws.config` | Production, AWS Batch queues |
| `test` | `conf/test.config` | CI/CD, TCGA-EE-A3J5 demo patient |

---

## Service Architecture (planned)

The goal is to wrap the pipeline in a **FastAPI service** so that:
1. Clinicians or bioinformaticians can submit a patient's data via API or web UI
2. The server queues and runs the pipeline asynchronously
3. Results are returned as structured JSON + a visual report
4. Multiple patients can be processed concurrently

### Target architecture

```
┌─────────────────────────────────────────────────────────────┐
│                        Frontend (React)                      │
│  Patient upload  │  Job status  │  Results viewer  │  Admin  │
└────────────────────────────┬────────────────────────────────┘
                             │ REST / WebSocket
                             ▼
┌─────────────────────────────────────────────────────────────┐
│                      FastAPI Server                          │
│                                                              │
│  POST /api/v1/jobs          Submit new patient job            │
│  GET  /api/v1/jobs/:id      Job status + progress            │
│  GET  /api/v1/jobs/:id/results   Download results            │
│  GET  /api/v1/jobs/:id/report    Vaccine summary report      │
│  POST /api/v1/predict       Lightweight: MAF → top 20 only   │
│                                                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐       │
│  │ Job Queue    │  │ Auth (JWT)   │  │ File Storage │       │
│  │ (Celery/RQ)  │  │              │  │ (S3)         │       │
│  └──────┬───────┘  └──────────────┘  └──────────────┘       │
│         │                                                    │
└─────────┼────────────────────────────────────────────────────┘
          │
          ▼
┌─────────────────────────────────────────────────────────────┐
│                    Pipeline Workers                           │
│                                                              │
│  Option A: Nextflow subprocess                               │
│    FastAPI spawns `nextflow run main.nf` per job             │
│    Pros: reuses all existing pipeline code                   │
│    Cons: heavyweight, Nextflow JVM startup, hard to stream   │
│                                                              │
│  Option B: Native Python pipeline                            │
│    Rewrite bin/ scripts as importable modules                │
│    Celery workers execute stages directly                    │
│    Pros: lightweight, fine-grained progress, no JVM          │
│    Cons: loses Nextflow resume/retry, need own orchestration │
│                                                              │
│  Option C: Hybrid (recommended)                              │
│    "Fast path" (modules 8-11): native Python, ~2 min         │
│    "Full path" (modules 1-11): Nextflow on AWS Batch         │
│    Pros: fast for common case, full pipeline when needed     │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Recommended approach: Hybrid (Option C)

Most users will submit a MAF + expression + HLA alleles (the output of standard bioinformatics pipelines like nf-core/sarek). For this case, we only need modules 8-11, which run in ~2 minutes and don't need Nextflow overhead.

**Fast path (native Python, modules 8-11):**
```python
# All bin/ scripts are already importable
from maf_to_pipeline_input import generate_candidates
from score_immunogenicity import score_candidates
from rank_and_select import rank_candidates
from design_polyepitope import design_polyepitope
from design_mrna import design_mrna
```

**Full path (Nextflow, modules 1-11):**
- For users submitting FASTQ/BAM, kick off Nextflow on AWS Batch
- Poll for completion, stream logs via WebSocket
- Return results when done

### API design

```
POST /api/v1/jobs
  Body: {
    patient_id: "PT-001",
    maf_file: <upload>,              # or S3 URI
    expression_file: <upload>,
    hla_alleles: ["HLA-A*02:01", "HLA-A*03:01", ...],
    tumor_type: "SKCM",
    profile: "high_tmb",             # or "low_tmb", "research"
    top_n: 20
  }
  Response: { job_id: "abc123", status: "queued" }

GET /api/v1/jobs/abc123
  Response: {
    status: "running",
    stage: "scoring",                # current pipeline stage
    progress: 0.6,
    started_at: "2026-04-07T12:00:00Z"
  }

GET /api/v1/jobs/abc123/results
  Response: {
    selected_neoantigens: [...],     # top 20 with scores
    polyepitope_sequence: "METDTL...",
    mrna_sequence: "AUGGAG...",
    synthesis_spec: {...},
    report_url: "/api/v1/jobs/abc123/report"
  }

POST /api/v1/predict
  # Lightweight endpoint: skip mRNA design, just return ranked neoantigens
  Body: { maf: <upload>, expression: <upload>, hla_alleles: [...] }
  Response: { candidates: [...top 20 with scores...] }
```

### Frontend (future)

A React/Next.js dashboard with:

1. **Upload page**: drag-and-drop MAF + expression + HLA, select tumor type and profile
2. **Job queue**: list of running/completed jobs with real-time progress
3. **Results viewer**:
   - Ranked neoantigen table (sortable, filterable)
   - Filter funnel visualization (Sankey diagram)
   - Per-epitope detail cards (binding, immunogenicity, structure, flags)
   - Polyepitope construct diagram with linkers highlighted
   - mRNA stats (length, GC%, MFE, CAI)
4. **Vaccine summary card**: one-page printable PDF per patient
5. **Admin**: user management, job history, compute cost tracking

### Migration steps

Rough order of work to get from current state to a working service:

1. **Refactor bin/ for importability**: make each script's `main()` callable with Python args (most are already close — they use argparse but also have importable functions)
2. **FastAPI skeleton**: `/predict` endpoint that runs modules 8-11 in-process
3. **Job queue**: add Celery + Redis for async job execution
4. **File handling**: S3 upload/download for input files and results
5. **Auth**: JWT-based auth (start with API keys for internal use)
6. **Frontend**: React dashboard with upload + results viewer
7. **Full pipeline**: wire Nextflow execution for FASTQ/BAM inputs
8. **Deployment**: Docker Compose for dev, ECS/Fargate for production

### Key decisions still open

- **Database**: PostgreSQL for job metadata, or DynamoDB for serverless?
- **Worker scaling**: fixed GPU instance vs. spot fleet for BigMHC?
- **Multi-tenancy**: single-tenant (one deployment per lab) or multi-tenant?
- **Model serving**: serve BigMHC via TorchServe/Triton for batched GPU inference, or keep in-process?
- **HIPAA**: if handling real patient data, need encryption at rest, audit logging, BAA with AWS

---

## Data flow

```
Patient data (MAF/VCF + expression + HLA)
    │
    ▼ uploaded to S3 or sent via API
    │
FastAPI server validates input
    │
    ▼ queues job
    │
Celery worker picks up job
    │
    ├── Fast path: runs modules 8-11 in Python (~2 min)
    │   └── writes results to S3
    │
    └── Full path: launches Nextflow on AWS Batch (~4-8 hr)
        └── polls for completion, writes results to S3
    │
    ▼
Results available via API / frontend
    │
    ▼
Vaccine summary card (HTML/PDF) generated on demand
```

---

## Current AWS resources

| Service | Resource | Details |
|---|---|---|
| EC2 | chatnav-gpu (g6e.4xlarge) | Benchmarking + BigMHC scoring, us-east-1 |
| EC2 | g6-big2 (g6.4xlarge) | General GPU, us-east-1 |
| EC2 | claude-cowork (t3.large) | Windows dev instance |
| S3 | s3://chatnav-pipeline-data/ | Pipeline data cache (~620GB) |
| S3 | s3://neoantigen-pipeline-workdir/ | Nextflow work directory |
| Batch | neoantigen-cpu-standard | m5.4xlarge compute environment |
| Batch | neoantigen-cpu-light | m5.xlarge compute environment |
| Batch | neoantigen-gpu-parabricks | g4dn.xlarge compute environment |
