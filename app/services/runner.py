from __future__ import annotations

import json
import logging
import os
import signal
import subprocess
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy.orm import Session

from app.config import Settings
from app.models import Job

logger = logging.getLogger("chatnav.runner")

# Keys that are API-only and should NOT be passed as nextflow --param flags
_SKIP_NEXTFLOW_KEYS = {"patient_id", "maf"}

# Valid nextflow params (from nextflow.config) — only these get passed as --flags
_NEXTFLOW_PARAMS = {
    "tumor_type", "tumor_fastq_1", "tumor_fastq_2",
    "normal_fastq_1", "normal_fastq_2", "tumor_bam", "normal_bam",
    "vcf", "cnv_segments", "fusions", "rnaseq_fastq_1", "rnaseq_fastq_2",
    "expression", "hla_alleles", "candidates_fasta", "candidates_meta",
    "binding_predictions", "parabricks", "structural_scoring",
    "weight_profile", "top_n", "min_ccf", "min_tpm",
    "binding_rank_threshold", "truncal_ccf",
}


def preprocess_maf(job: Job, settings: Settings) -> dict:
    """Run maf_to_pipeline_input.py to convert MAF → binding_predictions + candidates_meta.

    Returns a dict of extra nextflow params to merge into the job params.
    """
    params = json.loads(job.params_json)
    maf_path = params.get("maf")
    if not maf_path:
        return {}

    upload_dir = Path(maf_path).parent
    output_binding = str(upload_dir / "binding_predictions.tsv")
    output_meta = str(upload_dir / "candidates_meta.tsv")

    script = settings.pipeline_root / "bin" / "maf_to_pipeline_input.py"
    cmd = [
        sys.executable, str(script),
        "--maf", maf_path,
        "--output-binding", output_binding,
        "--output-meta", output_meta,
    ]

    # Pass optional args if available
    if params.get("expression"):
        cmd.extend(["--expression", params["expression"]])
    else:
        # expression is required by the script — create a minimal placeholder
        raise RuntimeError(
            "MAF submission requires an expression file. "
            "Upload an expression TSV or provide tumor_type for GTEx fallback."
        )

    if params.get("hla_alleles"):
        cmd.extend(["--hla-alleles", params["hla_alleles"]])
    else:
        raise RuntimeError(
            "MAF submission requires HLA alleles. "
            "Upload an hla_alleles file."
        )

    logger.info("Pre-processing MAF: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(settings.pipeline_root),
        timeout=3600,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"MAF pre-processing failed (exit {result.returncode}):\n"
            + result.stderr[-2000:]
        )

    return {
        "binding_predictions": output_binding,
        "candidates_meta": output_meta,
    }


def launch_nextflow(
    job: Job,
    settings: Settings,
) -> subprocess.Popen:
    """Build the nextflow command and launch it as a background subprocess."""
    params = json.loads(job.params_json)

    cmd = [
        settings.nextflow_bin,
        "run",
        str(settings.pipeline_root / "main.nf"),
        "-profile",
        settings.nextflow_profile,
        "--patient_id",
        job.patient_id,
        "--outdir",
        str(settings.data_dir / "results" / job.id),
    ]

    # Only pass recognized nextflow params
    for key, value in params.items():
        if key in _SKIP_NEXTFLOW_KEYS or key not in _NEXTFLOW_PARAMS:
            continue
        if value is None:
            continue
        cmd.extend([f"--{key}", str(value)])

    log_dir = settings.data_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{job.id}.log"
    log_file = open(log_path, "w")  # noqa: SIM115

    logger.info("Launching nextflow: %s", " ".join(cmd))

    proc = subprocess.Popen(
        cmd,
        stdout=log_file,
        stderr=subprocess.STDOUT,
        cwd=str(settings.pipeline_root),
        start_new_session=True,
    )
    return proc


def _read_tail(path: Path, lines: int = 50) -> str:
    """Read last N lines of a file, returning empty string if missing."""
    try:
        text = path.read_text()
        return "\n".join(text.splitlines()[-lines:])
    except OSError:
        return ""


class JobMonitor:
    """Background thread that polls running Nextflow processes for completion."""

    def __init__(self, settings: Settings, session_factory):
        self._active_jobs: dict[str, subprocess.Popen] = {}
        self._stop_event = threading.Event()
        self._settings = settings
        self._session_factory = session_factory
        self._lock = threading.Lock()
        self._thread: threading.Thread | None = None

    def register(self, job_id: str, proc: subprocess.Popen):
        """Register a newly launched process for monitoring."""
        with self._lock:
            self._active_jobs[job_id] = proc

    def recover_orphaned_jobs(self):
        """On startup, mark any 'running' jobs with dead PIDs as failed."""
        with self._session_factory() as session:
            running = (
                session.query(Job).filter(Job.status == "running").all()
            )
            for job in running:
                if job.pid and self._pid_alive(job.pid):
                    logger.info(
                        "Job %s (PID %d) still alive, will not monitor subprocess directly",
                        job.id,
                        job.pid,
                    )
                else:
                    logger.warning(
                        "Job %s has no live process — marking failed", job.id
                    )
                    job.status = "failed"
                    job.error_message = (
                        "Process not found on startup (API restart or crash)"
                    )
                    job.completed_at = datetime.now(timezone.utc)
            session.commit()

    def _launch_pending(self):
        """Launch pending jobs if there are free slots."""
        with self._lock:
            active_count = len(self._active_jobs)
        if active_count >= self._settings.max_concurrent_jobs:
            return

        with self._session_factory() as session:
            pending = (
                session.query(Job)
                .filter(Job.status == "pending")
                .order_by(Job.created_at)
                .limit(self._settings.max_concurrent_jobs - active_count)
                .all()
            )
            for job in pending:
                try:
                    # If MAF was uploaded, pre-process it first
                    params = json.loads(job.params_json)
                    if params.get("maf"):
                        extra = preprocess_maf(job, self._settings)
                        params.update(extra)
                        job.params_json = json.dumps(params)

                    proc = launch_nextflow(job, self._settings)
                    job.status = "running"
                    job.started_at = datetime.now(timezone.utc)
                    job.pid = proc.pid
                    self.register(job.id, proc)
                    logger.info("Started pending job %s (PID %d)", job.id, proc.pid)
                except Exception as e:
                    logger.exception("Failed to launch job %s", job.id)
                    job.status = "failed"
                    job.error_message = str(e)
                    job.completed_at = datetime.now(timezone.utc)
            session.commit()

    def _poll_loop(self):
        while not self._stop_event.is_set():
            # Check active jobs for completion
            with self._lock:
                jobs_snapshot = list(self._active_jobs.items())

            for job_id, proc in jobs_snapshot:
                retcode = proc.poll()
                if retcode is None:
                    continue

                log_path = self._settings.data_dir / "logs" / f"{job_id}.log"
                with self._session_factory() as session:
                    job = session.get(Job, job_id)
                    if job is None:
                        continue
                    if retcode == 0:
                        job.status = "completed"
                        logger.info("Job %s completed successfully", job_id)
                    else:
                        job.status = "failed"
                        job.error_message = (
                            f"Nextflow exited with code {retcode}\n"
                            + _read_tail(log_path)
                        )
                        logger.warning("Job %s failed (exit %d)", job_id, retcode)
                    job.completed_at = datetime.now(timezone.utc)
                    session.commit()

                with self._lock:
                    self._active_jobs.pop(job_id, None)

            # Launch any pending jobs that are waiting for a slot
            self._launch_pending()

            self._stop_event.wait(5.0)

    @staticmethod
    def _pid_alive(pid: int) -> bool:
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False

    def cancel_job(self, job_id: str) -> bool:
        """Send SIGTERM to a running job. Returns True if signal was sent."""
        with self._lock:
            proc = self._active_jobs.get(job_id)

        if proc is not None:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                return True
            except OSError:
                return False

        # Fallback: try PID from DB
        with self._session_factory() as session:
            job = session.get(Job, job_id)
            if job and job.pid:
                try:
                    os.killpg(os.getpgid(job.pid), signal.SIGTERM)
                    return True
                except OSError:
                    return False
        return False

    def start(self):
        self._thread = threading.Thread(target=self._poll_loop, daemon=True)
        self._thread.start()
        logger.info("JobMonitor started")

    def stop(self):
        self._stop_event.set()
        if self._thread:
            self._thread.join(timeout=10)
        logger.info("JobMonitor stopped")
