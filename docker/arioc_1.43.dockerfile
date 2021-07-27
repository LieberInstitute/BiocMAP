#  This likely needs a lot of work! For example, making sure 'unzip', CUDA,
#  and other dependencies exist. Are we assuming 'create_images.sh' has access
#  to an Nvidia GPU?

FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

RUN mkdir arioc && \
    cd arioc && \
    wget https://github.com/RWilton/Arioc/releases/download/v1.43/Arioc.x.143.zip && \
    unzip Arioc.x.143.zip
        
RUN cd src && \
    make clean && \
    make AriocE && \
    make AriocU && \
    make AriocP && \
    cp ../bin/* /usr/local/bin/
