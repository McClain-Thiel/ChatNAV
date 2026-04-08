from __future__ import annotations

import csv
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path

from fastapi import (
    APIRouter,
    Depends,
    File,
    Form,
    HTTPException,
    Query,
    Request,
    UploadFile,
)
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from app.config import settings
from app.database import get_db
from app.models import Job
from app.schemas import (
    JobDetailResponse,
    JobListResponse,
    JobOutputsResponse,
    JobResponse,
    LogResponse,
)
from app.services.files import (
    guess_content_type,
    list_outputs,
    resolve_output_path,
    save_upload,
)
from app.services.pipeline import (
    build_modules_state,
    build_vaccine_result,
)
from app.services.runner import launch_nextflow

router = APIRouter(prefix="/api/v1/jobs", tags=["jobs"])


def _job_to_response(job: Job) -> JobResponse:
    return JobResponse(
        id=job.id,
        patient_id=job.patient_id,
        status=job.status,
        entry_point=job.entry_point,
        created_at=job.created_at,
        started_at=job.started_at,
        completed_at=job.completed_at,
        error_message=job.error_message,
        params=json.loads(job.params_json),
    )


def _get_job_or_404(job_id: str, db: Session) -> Job:
    job = db.get(Job, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


# ── File upload param names that map to nextflow --param flags ──
FILE_PARAMS = [
    "vcf",
    "maf",
    "hla_alleles",
    "expression",
    "cnv_segments",
    "tumor_fastq_1",
    "tumor_fastq_2",
    "normal_fastq_1",
    "normal_fastq_2",
    "rnaseq_fastq_1",
    "rnaseq_fastq_2",
    "candidates_fasta",
    "candidates_meta",
    "binding_predictions",
]


@router.post("", status_code=202)
async def submit_job(
    request: Request,
    patient_id: str = Form(...),
    tumor_type: str | None = Form(None),
    weight_profile: str = Form("high_tmb"),
    top_n: int = Form(20),
    structural_scoring: bool = Form(False),
    vcf: UploadFile | None = File(None),
    maf: UploadFile | None = File(None),
    hla_alleles: UploadFile | None = File(None),
    expression: UploadFile | None = File(None),
    cnv_segments: UploadFile | None = File(None),
    tumor_fastq_1: UploadFile | None = File(None),
    tumor_fastq_2: UploadFile | None = File(None),
    normal_fastq_1: UploadFile | None = File(None),
    normal_fastq_2: UploadFile | None = File(None),
    rnaseq_fastq_1: UploadFile | None = File(None),
    rnaseq_fastq_2: UploadFile | None = File(None),
    candidates_fasta: UploadFile | None = File(None),
    candidates_meta: UploadFile | None = File(None),
    binding_predictions: UploadFile | None = File(None),
    db: Session = Depends(get_db),
) -> JobResponse:
    job_id = str(uuid.uuid4())
    upload_dir = settings.data_dir / "uploads" / job_id

    # Build nextflow params dict
    params: dict = {
        "weight_profile": weight_profile,
        "top_n": top_n,
    }
    if tumor_type:
        params["tumor_type"] = tumor_type
    if structural_scoring:
        params["structural_scoring"] = "true"

    # Save uploaded files and add their paths to params
    file_uploads = {
        "vcf": vcf,
        "maf": maf,
        "hla_alleles": hla_alleles,
        "expression": expression,
        "cnv_segments": cnv_segments,
        "tumor_fastq_1": tumor_fastq_1,
        "tumor_fastq_2": tumor_fastq_2,
        "normal_fastq_1": normal_fastq_1,
        "normal_fastq_2": normal_fastq_2,
        "rnaseq_fastq_1": rnaseq_fastq_1,
        "rnaseq_fastq_2": rnaseq_fastq_2,
        "candidates_fasta": candidates_fasta,
        "candidates_meta": candidates_meta,
        "binding_predictions": binding_predictions,
    }

    for param_name, upload_file in file_uploads.items():
        if upload_file is not None and upload_file.filename:
            saved_path = await save_upload(upload_file, upload_dir, param_name)
            params[param_name] = str(saved_path)

    # Create job record
    job = Job(
        id=job_id,
        patient_id=patient_id,
        status="pending",
        params_json=json.dumps(params),
        created_at=datetime.now(timezone.utc),
    )
    db.add(job)
    db.commit()
    db.refresh(job)

    # The JobMonitor background thread will pick up pending jobs and launch them
    return _job_to_response(job)


@router.get("")
def list_jobs(
    status: str | None = Query(None),
    limit: int = Query(20, ge=1, le=100),
    offset: int = Query(0, ge=0),
    db: Session = Depends(get_db),
) -> JobListResponse:
    query = db.query(Job)
    if status:
        query = query.filter(Job.status == status)
    total = query.count()
    jobs = query.order_by(Job.created_at.desc()).offset(offset).limit(limit).all()
    return JobListResponse(
        jobs=[_job_to_response(j) for j in jobs],
        total=total,
    )


@router.get("/{job_id}")
def get_job(
    job_id: str,
    db: Session = Depends(get_db),
) -> JobDetailResponse:
    job = _get_job_or_404(job_id, db)
    results_dir = settings.data_dir / "results" / job.id
    modules = build_modules_state(job.status, results_dir)

    # Determine current module (first non-complete, or last)
    current_module = 1
    for m in modules:
        if m["status"] != "complete":
            current_module = m["number"]
            break
    else:
        current_module = modules[-1]["number"] if modules else 1

    # Build vaccine result if job is complete
    result = None
    if job.status == "completed":
        result = build_vaccine_result(results_dir, job.patient_id)

    return JobDetailResponse(
        id=job.id,
        patient_id=job.patient_id,
        status=job.status,
        entry_point=job.entry_point,
        created_at=job.created_at,
        started_at=job.started_at,
        completed_at=job.completed_at,
        error_message=job.error_message,
        params=json.loads(job.params_json),
        modules=modules,
        current_module=current_module,
        result=result,
    )


@router.get("/{job_id}/outputs")
def get_job_outputs(
    job_id: str,
    db: Session = Depends(get_db),
) -> JobOutputsResponse:
    job = _get_job_or_404(job_id, db)
    files = list_outputs(job_id, settings)
    return JobOutputsResponse(
        job_id=job.id,
        status=job.status,
        files=files,
    )


@router.get("/{job_id}/outputs/{file_path:path}")
def download_output(
    job_id: str,
    file_path: str,
    db: Session = Depends(get_db),
):
    _get_job_or_404(job_id, db)
    resolved = resolve_output_path(job_id, file_path, settings)
    return FileResponse(
        path=str(resolved),
        media_type=guess_content_type(resolved),
        filename=resolved.name,
    )


@router.get("/{job_id}/preview/{file_path:path}")
def preview_file(
    job_id: str,
    file_path: str,
    max_rows: int = Query(100, ge=1, le=5000),
    db: Session = Depends(get_db),
):
    """Return file content parsed for frontend display.

    TSV/CSV → {"type": "table", "columns": [...], "rows": [[...], ...]}
    JSON    → {"type": "json", "data": {...}}
    FASTA   → {"type": "fasta", "sequences": [{"header": "...", "sequence": "..."}]}
    Other   → {"type": "text", "content": "..."}
    """
    _get_job_or_404(job_id, db)
    resolved = resolve_output_path(job_id, file_path, settings)
    name = resolved.name.lower()

    if name.endswith(".tsv") or name.endswith(".csv"):
        delimiter = "\t" if name.endswith(".tsv") else ","
        rows = []
        columns = []
        with open(resolved, newline="") as f:
            reader = csv.reader(f, delimiter=delimiter)
            for i, row in enumerate(reader):
                if i == 0:
                    columns = row
                elif i <= max_rows:
                    rows.append(row)
        return {"type": "table", "columns": columns, "rows": rows, "total_rows": len(rows)}

    if name.endswith(".json"):
        data = json.loads(resolved.read_text())
        return {"type": "json", "data": data}

    if name.endswith(".fasta") or name.endswith(".faa") or name.endswith(".fa"):
        text = resolved.read_text()
        sequences = []
        header = ""
        seq_lines: list[str] = []
        for line in text.splitlines():
            if line.startswith(">"):
                if header or seq_lines:
                    sequences.append({"header": header, "sequence": "".join(seq_lines)})
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header or seq_lines:
            sequences.append({"header": header, "sequence": "".join(seq_lines)})
        return {"type": "fasta", "sequences": sequences}

    # Plain text fallback
    text = resolved.read_text()
    if len(text) > 50_000:
        text = text[:50_000] + "\n... (truncated)"
    return {"type": "text", "content": text}


@router.post("/{job_id}/cancel")
def cancel_job(
    job_id: str,
    request: Request,
    db: Session = Depends(get_db),
) -> JobResponse:
    job = _get_job_or_404(job_id, db)
    if job.status not in ("pending", "running"):
        raise HTTPException(
            status_code=409,
            detail=f"Cannot cancel job with status '{job.status}'",
        )

    monitor = request.app.state.monitor
    if job.status == "running":
        monitor.cancel_job(job_id)

    job.status = "cancelled"
    job.completed_at = datetime.now(timezone.utc)
    db.commit()
    db.refresh(job)
    return _job_to_response(job)


@router.get("/{job_id}/log")
def get_job_log(
    job_id: str,
    tail: int = Query(0, ge=0, description="Return last N lines (0 = all)"),
    db: Session = Depends(get_db),
) -> LogResponse:
    _get_job_or_404(job_id, db)
    log_path = settings.data_dir / "logs" / f"{job_id}.log"
    if not log_path.is_file():
        return LogResponse(job_id=job_id, log="")

    text = log_path.read_text()
    if tail > 0:
        text = "\n".join(text.splitlines()[-tail:])
    return LogResponse(job_id=job_id, log=text)
