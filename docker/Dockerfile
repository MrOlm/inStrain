# Base Image
FROM continuumio/miniconda3:4.6.14
USER root:root

# Not sure why, but errors out without
USER root:root
RUN mkdir -p /mnt
WORKDIR /mnt
#RUN apt-get --allow-releaseinfo-change update

# Add channels
RUN conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge

#RUN conda install mamba
RUN conda create -n work -y instrain awscli samtools coverm && conda clean -afy

RUN /bin/bash -c "source activate work && \
    pip install boto3 && pip install instrain --upgrade"

COPY . .
RUN chmod +rx *

# Metadata
LABEL container.maintainer="Matt Olm <mattolm@stanford.edu>" \
      container.base.image="continuumio/miniconda3" \ 
      software.name="inStrain"
