FROM python:3.11-slim

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 6: Candidate generation (pVACtools)"
LABEL version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    pvactools>=4.0 \
    pandas>=2.0

WORKDIR /work
ENTRYPOINT []
