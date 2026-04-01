FROM ubuntu:22.04

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 4: Expression quantification (Salmon)"
LABEL version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        python3 \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Salmon 1.10
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz \
    && tar -xzf salmon-1.10.0_linux_x86_64.tar.gz -C /opt/ \
    && rm salmon-1.10.0_linux_x86_64.tar.gz
ENV PATH="/opt/salmon-latest_linux_x86_64/bin:${PATH}"

RUN pip3 install --no-cache-dir pandas

WORKDIR /work
ENTRYPOINT []
