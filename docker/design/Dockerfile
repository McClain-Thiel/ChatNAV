FROM python:3.11-slim

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 10-11: Polyepitope + mRNA design"
LABEL version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    pandas>=2.0 \
    numpy>=1.24 \
    biopython>=1.81

# ── LinearDesign (Baidu, 2023) ──
# Joint codon + mRNA structure optimization
# Clone and compile: https://github.com/LinearDesignSoftware/LinearDesign
# RUN git clone https://github.com/LinearDesignSoftware/LinearDesign.git /opt/LinearDesign \
#     && cd /opt/LinearDesign \
#     && make
# ENV PATH="/opt/LinearDesign:${PATH}"
#
# If LinearDesign is not available, the pipeline falls back to CAI
# (Codon Adaptation Index) maximization using human codon usage tables.

WORKDIR /work

ENTRYPOINT []
