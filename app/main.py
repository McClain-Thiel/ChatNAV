from contextlib import asynccontextmanager
import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.config import settings
from app.database import SessionLocal, create_db_tables
from app.routes.jobs import router as jobs_router
from app.services.demo import seed_demo_job
from app.services.runner import JobMonitor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
)
logger = logging.getLogger("chatnav")


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    settings.data_dir.mkdir(parents=True, exist_ok=True)
    create_db_tables()

    seed_demo_job(settings, SessionLocal)

    monitor = JobMonitor(settings, SessionLocal)
    monitor.recover_orphaned_jobs()
    monitor.start()
    app.state.monitor = monitor

    logger.info("ChatNAV API started — pipeline root: %s", settings.pipeline_root)
    yield

    # Shutdown
    monitor.stop()
    logger.info("ChatNAV API stopped")


app = FastAPI(
    title="ChatNAV",
    description="Personalized mRNA neoantigen vaccine design API",
    version="0.1.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(jobs_router)
