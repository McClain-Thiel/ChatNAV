FROM ubuntu:22.04

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Modules 1-2: Preprocessing (BWA-MEM2, GATK, samtools)"
LABEL version="0.1.0"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        curl \
        unzip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        openjdk-17-jre-headless \
        python3 \
        tabix \
        bcftools \
    && rm -rf /var/lib/apt/lists/*

# BWA-MEM2 v2.2.1
RUN wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    && tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2 -C /opt/ \
    && rm bwa-mem2-2.2.1_x64-linux.tar.bz2
ENV PATH="/opt/bwa-mem2-2.2.1_x64-linux:${PATH}"

# samtools 1.19
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 \
    && tar -xjf samtools-1.19.tar.bz2 \
    && cd samtools-1.19 && ./configure && make -j$(nproc) && make install \
    && cd .. && rm -rf samtools-1.19*

# GATK 4.5.0.0
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
    && unzip -q gatk-4.5.0.0.zip -d /opt/ \
    && rm gatk-4.5.0.0.zip
ENV PATH="/opt/gatk-4.5.0.0:${PATH}"

WORKDIR /work
ENTRYPOINT []
