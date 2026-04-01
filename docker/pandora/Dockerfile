FROM continuumio/miniconda3:latest

RUN conda install -y -c salilab -c conda-forge modeller && \
    pip install CSB-PANDORA biopython pdb2sql && \
    conda clean -afy

# Set MODELLER license key (must be done before import)
RUN sed -i "s/license = 'XXXX'/license = 'MODELIRANJE'/" \
    /opt/conda/lib/modeller-10.8/modlib/modeller/config.py && \
    python -c "import modeller; print('MODELLER OK')"

WORKDIR /work
