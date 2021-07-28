FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

##  hisat2 (2.2.1)
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.tar.gz && \
    tar -xzf v2.2.1.tar.gz && \
    cd hisat2-2.2.1 && \
    make && \
    cp hisat2 /usr/local/bin/ && \
    cp hisat2-align* /usr/local/bin/ && \
    cp hisat2-build* /usr/local/bin/ && \
    cp hisat2-inspect* /usr/local/bin/ && \
    cp *.py /usr/local/bin/ && \
    cd /usr/local/src
        
##  bismark itself (0.23.0) 
RUN wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz && \
    tar -xzf 0.23.0.tar.gz && \
    cp Bismark-0.23.0/bismark* /usr/local/bin/ && \
    cp Bismark-0.23.0/coverage2cytosine /usr/local/bin/
