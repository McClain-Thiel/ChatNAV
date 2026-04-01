FROM nvcr.io/nvidia/clara/clara-parabricks:4.1.0-1

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Modules 1-2: GPU-accelerated preprocessing (Parabricks)"
LABEL version="0.1.0"

# Parabricks is free for academic/research use (< $5M revenue).
# Commercial use requires NVIDIA license (~$10K+/yr).
# See: https://www.nvidia.com/en-us/clara/parabricks/

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
ENTRYPOINT []
