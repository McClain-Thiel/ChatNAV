FROM ubuntu:22.04

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 3: HLA typing (HLA-HD v1.7)"
LABEL version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        bowtie2 \
        samtools \
        python3 \
        python3-pip \
        wget \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pandas

# ══════════════════════════════════════════════════════════════════
# ACADEMIC LICENSE REQUIRED
#
# HLA-HD v1.7 requires a free academic license.
# Download from: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
# Place the tarball in this directory before building.
# ══════════════════════════════════════════════════════════════════
COPY hlahd.*.tar.gz /tmp/
RUN tar -xzf /tmp/hlahd.*.tar.gz -C /opt/ && \
    cd /opt/hlahd* && \
    sh install.sh && \
    rm /tmp/hlahd.*.tar.gz
ENV PATH="/opt/hlahd/bin:${PATH}"

WORKDIR /work
ENTRYPOINT []
