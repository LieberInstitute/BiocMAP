FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

#  Install samtools (1.10)
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjf samtools.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install && \
    cd /usr/local/src
        
##  bismark itself (0.23.0) 
RUN wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz && \
    tar -xzf 0.23.0.tar.gz && \
    cp Bismark-0.23.0/bismark* /usr/local/bin/ && \
    cp Bismark-0.23.0/coverage2cytosine /usr/local/bin/
