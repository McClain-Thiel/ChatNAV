FROM trinityctat/starfusion:1.12.0

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 2: Fusion detection (STAR-Fusion)"
LABEL version="0.1.0"

# STAR-Fusion genome library is NOT baked into the image (30GB+).
# Download from: https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
# Mount at runtime: -v /path/to/ctat_genome_lib:/genome_lib

RUN pip install --no-cache-dir pandas

WORKDIR /work
