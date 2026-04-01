FROM ensemblorg/ensembl-vep:release_110.1

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 2: Variant annotation (Ensembl VEP)"
LABEL version="0.1.0"

# VEP cache is NOT baked into the image (50GB+).
# Mount it at runtime: -v /path/to/vep_cache:/opt/vep/.vep
# Download with: vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38

USER root
RUN apt-get update && apt-get install -y --no-install-recommends tabix && rm -rf /var/lib/apt/lists/*
USER vep

WORKDIR /work
