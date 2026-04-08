from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_prefix="CHATNAV_")

    # Paths
    pipeline_root: Path = Path(__file__).resolve().parent.parent
    data_dir: Path = Path(__file__).resolve().parent.parent / "data"
    database_url: str = f"sqlite:///{Path(__file__).resolve().parent.parent / 'data' / 'chatnav.db'}"
    nextflow_bin: str = "nextflow"
    nextflow_profile: str = "local"

    # Limits
    max_upload_size_mb: int = 500
    max_concurrent_jobs: int = 3


settings = Settings()
