FROM nvidia/cuda:11.4.1-devel-ubuntu18.04

WORKDIR /usr/local/src

RUN apt-get update
RUN apt-get install -y wget make zip

RUN mkdir arioc && \
    cd arioc && \
    wget https://github.com/RWilton/Arioc/releases/download/v1.43/Arioc.x.143.zip && \
    unzip Arioc.x.143.zip
        
RUN cd arioc/src && \
    make clean && \
    make AriocE && \
    make AriocU && \
    make AriocP && \
    cp ../bin/* /usr/local/bin/
