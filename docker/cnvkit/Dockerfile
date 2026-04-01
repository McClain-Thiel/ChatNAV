FROM python:3.11-slim

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 2: CNV calling (CNVkit)"
LABEL version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base-core \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    cnvkit>=0.9.10 \
    pandas numpy

WORKDIR /work
ENTRYPOINT []
