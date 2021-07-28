FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

ENV LD_LIBRARY_PATH=/usr/local/lib

#  Install samtools (1.10)
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjf samtools.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install && \
    cd /usr/local/src
    
##  Install libBigWig (0.4.6), a MethylDackel dependency
RUN wget https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.6.tar.gz && \
    tar -xzf 0.4.6.tar.gz && \
    cd libBigWig-0.4.6 && \
    make install && \
    cd /usr/local/src
          
##  Install htslib (1.12), a MethylDackel dependency
RUN wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
    tar -xjf htslib-1.12.tar.bz2 && \
    cd htslib-1.12 && \
    ./configure && \
    make && \
    make install && \
    cd /usr/local/src
        
##  MethylDackel itself (0.5.2)
RUN wget https://github.com/dpryan79/MethylDackel/archive/refs/tags/0.5.2.tar.gz && \
    tar -xzf 0.5.2.tar.gz && \
    cd MethylDackel-0.5.2 && \
    make LIBBIGWIG="/usr/local/lib/libBigWig.a" CFLAGS="-O3 -Wall -I/usr/local/include " LIBS="-L/usr/local/lib" && \
    make LIBBIGWIG="/usr/local/lib/libBigWig.a" CFLAGS="-O3 -Wall -I/usr/local/include " LIBS="-L/usr/local/lib" install
