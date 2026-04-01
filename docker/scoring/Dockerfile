FROM python:3.11-slim

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 8-9: Immunogenicity scoring and ranking"
LABEL version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    pandas>=2.0 \
    numpy>=1.24 \
    biopython>=1.81 \
    scipy>=1.11 \
    pyyaml>=6.0 \
    requests>=2.31 \
    tensorflow>=2.13

# DeepImmuno-CNN — open source immunogenicity prediction
# https://github.com/frankligy/DeepImmuno
RUN git clone https://github.com/frankligy/DeepImmuno.git /opt/DeepImmuno
ENV DEEPIMMUNO_MODEL_DIR=/opt/DeepImmuno/models

WORKDIR /work
ENTRYPOINT []
