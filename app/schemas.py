from datetime import datetime

from pydantic import BaseModel


class JobResponse(BaseModel):
    id: str
    patient_id: str
    status: str
    entry_point: int | None = None
    created_at: datetime
    started_at: datetime | None = None
    completed_at: datetime | None = None
    error_message: str | None = None
    params: dict

    model_config = {"from_attributes": True}


class JobListResponse(BaseModel):
    jobs: list[JobResponse]
    total: int


class OutputFile(BaseModel):
    path: str
    size_bytes: int
    modified_at: datetime


class JobOutputsResponse(BaseModel):
    job_id: str
    status: str
    files: list[OutputFile]


class LogResponse(BaseModel):
    job_id: str
    log: str


# ── Enriched job detail for frontend ──

class JobDetailResponse(BaseModel):
    """Full job detail including pipeline module states and vaccine results."""
    id: str
    patient_id: str
    status: str
    entry_point: int | None = None
    created_at: datetime
    started_at: datetime | None = None
    completed_at: datetime | None = None
    error_message: str | None = None
    params: dict
    modules: list[dict]
    current_module: int
    result: dict | None = None
