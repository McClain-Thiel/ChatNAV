FROM python:3.11-slim

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 7: MHC binding (MHCflurry + MHCnuggets)"
LABEL version="0.1.0"
LABEL license="Apache-2.0 (MHCflurry) + MIT (MHCnuggets)"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# MHCflurry 2.0 — Class I binding + antigen processing (Apache 2.0)
# MHCnuggets — Class II binding (MIT)
RUN pip install --no-cache-dir \
    mhcflurry>=2.1.0 \
    mhcnuggets>=2.4.0 \
    pandas>=2.0 \
    numpy>=1.24

# Download MHCflurry models (~2GB, needed at runtime)
RUN mhcflurry-downloads fetch models_class1_presentation

WORKDIR /work
ENTRYPOINT []
