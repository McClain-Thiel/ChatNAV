"""Seed a demo job from existing pipeline results.

Creates a pre-completed job in the database pointing to the real TCGA_EE_A3J5
pipeline output so users can explore the full results view immediately.
"""

from __future__ import annotations

import json
import logging
import shutil
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy.orm import Session

from app.config import Settings
from app.models import Job

logger = logging.getLogger("chatnav.demo")

DEMO_JOB_ID = "demo"
DEMO_PATIENT_ID = "Demo — TCGA Melanoma"

# Source results from the main repo's completed run
DEMO_RESULTS_SOURCE = Path("/home/ubuntu/ChatNAV/results/TCGA_EE_A3J5")


def seed_demo_job(settings: Settings, session_factory) -> None:
    """Create the demo job if it doesn't already exist."""
    with session_factory() as session:
        existing = session.get(Job, DEMO_JOB_ID)
        if existing is not None:
            logger.info("Demo job already exists, skipping seed")
            return

        # Copy results into the expected output directory structure
        dest = settings.data_dir / "results" / DEMO_JOB_ID
        if not dest.exists() and DEMO_RESULTS_SOURCE.exists():
            dest.mkdir(parents=True, exist_ok=True)
            # The results are flat — copy them into module subdirectories
            # Pipeline stats
            _copy_if_exists(DEMO_RESULTS_SOURCE / "pipeline_stats.json",
                           dest / "pipeline_stats.json")

            # Module 8: immunogenicity_scores.tsv
            (dest / "module_08").mkdir(exist_ok=True)
            _copy_if_exists(DEMO_RESULTS_SOURCE / "immunogenicity_scores.tsv",
                           dest / "module_08" / "immunogenicity_scores.tsv")

            # Module 9: selected_neoantigens.tsv
            (dest / "module_09").mkdir(exist_ok=True)
            _copy_if_exists(DEMO_RESULTS_SOURCE / "selected_neoantigens.tsv",
                           dest / "module_09" / "selected_neoantigens.tsv")

            # Module 10: polyepitope.faa, polyepitope_design.tsv
            (dest / "module_10").mkdir(exist_ok=True)
            _copy_if_exists(DEMO_RESULTS_SOURCE / "polyepitope.faa",
                           dest / "module_10" / "polyepitope.faa")
            _copy_if_exists(DEMO_RESULTS_SOURCE / "polyepitope_design.tsv",
                           dest / "module_10" / "polyepitope_design.tsv")

            # Module 11: mrna_sequence.fasta, synthesis_spec.json
            (dest / "module_11").mkdir(exist_ok=True)
            _copy_if_exists(DEMO_RESULTS_SOURCE / "mrna_sequence.fasta",
                           dest / "module_11" / "mrna_sequence.fasta")
            _copy_if_exists(DEMO_RESULTS_SOURCE / "synthesis_spec.json",
                           dest / "module_11" / "synthesis_spec.json")

            logger.info("Copied demo results to %s", dest)
        elif not DEMO_RESULTS_SOURCE.exists():
            logger.warning("Demo results source not found: %s", DEMO_RESULTS_SOURCE)
            return

        job = Job(
            id=DEMO_JOB_ID,
            patient_id=DEMO_PATIENT_ID,
            status="completed",
            params_json=json.dumps({
                "tumor_type": "Melanoma (SKCM)",
                "weight_profile": "high_tmb",
                "top_n": 20,
            }),
            created_at=datetime(2026, 3, 30, 14, 22, 0, tzinfo=timezone.utc),
            started_at=datetime(2026, 3, 30, 14, 22, 5, tzinfo=timezone.utc),
            completed_at=datetime(2026, 3, 30, 15, 48, 30, tzinfo=timezone.utc),
        )
        session.add(job)
        session.commit()
        logger.info("Seeded demo job: %s (%s)", DEMO_JOB_ID, DEMO_PATIENT_ID)


def _copy_if_exists(src: Path, dst: Path) -> None:
    if src.is_file():
        shutil.copy2(src, dst)
