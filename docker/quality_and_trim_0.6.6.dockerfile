FROM ubuntu:20.04

#  Bypass user input prompts during installation of java
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /usr/local/src

RUN apt-get update
RUN apt-get install -y wget make zip perl default-jre python3-dev python3-pip python3-setuptools

#  Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    chmod -R 775 FastQC && \
    ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc
    
#  Install cutadapt
RUN pip install cutadapt

#  Install Trim Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz && \
    tar -xzf 0.6.6.tar.gz && \
    cp TrimGalore-0.6.6/trim_galore /usr/local/bin/
