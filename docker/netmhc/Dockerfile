FROM ubuntu:22.04

LABEL maintainer="McClain Thiel"
LABEL description="Neoantigen pipeline Module 7: MHC binding prediction suite"
LABEL version="0.1.0"

# NetMHCpan suite requires tcsh
RUN apt-get update && apt-get install -y --no-install-recommends \
        tcsh \
        gawk \
        wget \
        python3 \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pandas numpy

# ══════════════════════════════════════════════════════════════════
# ACADEMIC LICENSE REQUIRED
#
# All NetMHC tools require free academic licenses from DTU Health Tech.
# Download the Linux tarballs and place them in this directory before
# building:
#
# NetMHCpan 4.1:    https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
# NetMHCIIpan 4.3:  https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/
# NetMHCstabpan 1.0: https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/
# NetChop 3.1:      https://services.healthtech.dtu.dk/services/NetChop-3.1/
# NetCTLpan 1.1:    https://services.healthtech.dtu.dk/services/NetCTLpan-1.1/
# ══════════════════════════════════════════════════════════════════

# NetMHCpan 4.1
COPY netMHCpan-4.1*.tar.gz /tmp/
RUN mkdir -p /opt/netMHCpan-4.1 && \
    tar -xzf /tmp/netMHCpan-4.1*.tar.gz -C /opt/ && \
    rm /tmp/netMHCpan-4.1*.tar.gz
ENV NMHOME=/opt/netMHCpan-4.1
ENV PATH="/opt/netMHCpan-4.1:${PATH}"

# NetMHCIIpan 4.3
COPY netMHCIIpan-4.3*.tar.gz /tmp/
RUN mkdir -p /opt/netMHCIIpan-4.3 && \
    tar -xzf /tmp/netMHCIIpan-4.3*.tar.gz -C /opt/ && \
    rm /tmp/netMHCIIpan-4.3*.tar.gz
ENV PATH="/opt/netMHCIIpan-4.3:${PATH}"

# NetMHCstabpan 1.0
COPY netMHCstabpan-1.0*.tar.gz /tmp/
RUN mkdir -p /opt/netMHCstabpan-1.0 && \
    tar -xzf /tmp/netMHCstabpan-1.0*.tar.gz -C /opt/ && \
    rm /tmp/netMHCstabpan-1.0*.tar.gz
ENV PATH="/opt/netMHCstabpan-1.0:${PATH}"

# NetChop 3.1
COPY netchop-3.1*.tar.gz /tmp/
RUN mkdir -p /opt/netchop-3.1 && \
    tar -xzf /tmp/netchop-3.1*.tar.gz -C /opt/ && \
    rm /tmp/netchop-3.1*.tar.gz
ENV PATH="/opt/netchop-3.1:${PATH}"

# NetCTLpan 1.1
COPY netCTLpan-1.1*.tar.gz /tmp/
RUN mkdir -p /opt/netCTLpan-1.1 && \
    tar -xzf /tmp/netCTLpan-1.1*.tar.gz -C /opt/ && \
    rm /tmp/netCTLpan-1.1*.tar.gz
ENV PATH="/opt/netCTLpan-1.1:${PATH}"

# Set tmp directory for NetMHC tools
RUN mkdir -p /tmp/netmhc
ENV TMPDIR=/tmp/netmhc

WORKDIR /work
ENTRYPOINT []
