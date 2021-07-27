FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

#  Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    chmod -R 775 FastQC && \
    cp FastQC/fastqc /usr/local/bin/

#  Install cutadapt
RUN pip install cutadapt

#  Install Trim Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz && \
    tar -xzf 0.6.6.tar.gz && \
    cp TrimGalore-0.6.6/trim_galore /usr/local/bin/
