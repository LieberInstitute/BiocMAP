FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

#  Install samblaster (v.0.1.26)
RUN wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz && \
    tar -xzf samblaster-v.0.1.26.tar.gz && \
    cd samblaster-v.0.1.26 && \
    make && \
    cp samblaster /usr/local/bin/

#  Install samtools (1.10)
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjf samtools.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install
