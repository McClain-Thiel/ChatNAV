from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import aiofiles
from fastapi import HTTPException, UploadFile

from app.config import Settings
from app.schemas import OutputFile

# Map file extensions to MIME types for common pipeline outputs
CONTENT_TYPES = {
    ".tsv": "text/tab-separated-values",
    ".csv": "text/csv",
    ".json": "application/json",
    ".fasta": "text/plain",
    ".faa": "text/plain",
    ".fa": "text/plain",
    ".txt": "text/plain",
    ".html": "text/html",
    ".vcf": "text/plain",
    ".vcf.gz": "application/gzip",
    ".gz": "application/gzip",
    ".tbi": "application/octet-stream",
    ".bam": "application/octet-stream",
    ".bai": "application/octet-stream",
    ".yaml": "text/yaml",
    ".yml": "text/yaml",
}


async def save_upload(
    file: UploadFile, job_dir: Path, param_name: str
) -> Path:
    """Save an uploaded file to the job's upload directory. Returns saved path."""
    job_dir.mkdir(parents=True, exist_ok=True)
    dest = job_dir / f"{param_name}_{file.filename}"
    async with aiofiles.open(dest, "wb") as f:
        while chunk := await file.read(8192):
            await f.write(chunk)
    return dest


def list_outputs(job_id: str, settings: Settings) -> list[OutputFile]:
    """Walk the results directory and return metadata for each output file."""
    results_dir = settings.data_dir / "results" / job_id
    if not results_dir.exists():
        return []

    files = []
    for path in sorted(results_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(results_dir)
        stat = path.stat()
        files.append(
            OutputFile(
                path=str(rel),
                size_bytes=stat.st_size,
                modified_at=datetime.fromtimestamp(
                    stat.st_mtime, tz=timezone.utc
                ),
            )
        )
    return files


def resolve_output_path(
    job_id: str, file_path: str, settings: Settings
) -> Path:
    """Resolve and validate that file_path stays within the job's results dir.

    Raises HTTPException on traversal attempts or missing files.
    """
    base = (settings.data_dir / "results" / job_id).resolve()
    target = (base / file_path).resolve()
    if not str(target).startswith(str(base)):
        raise HTTPException(status_code=400, detail="Invalid file path")
    if not target.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    return target


def guess_content_type(path: Path) -> str:
    """Guess MIME type from file extension."""
    name = path.name.lower()
    # Check compound extensions first (e.g. .vcf.gz)
    for ext, ct in CONTENT_TYPES.items():
        if name.endswith(ext):
            return ct
    return "application/octet-stream"
