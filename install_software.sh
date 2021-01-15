#!/bin/bash

#  IMPORTANT: run this script inside the repository directory to ensure 
#             software locations are properly linked. See
#             conf/command_path_long.config if you are interested in manually
#             configuring different software paths.
#
#  Usage:  bash install_software.sh [installation_type]

set -e

#  Verify java can be executed, since this is a pre-requisite
if [ -x "$(command -v java)" ]; then
    echo "Found a java runtime. Proceeding with the setup..."
    
    INSTALL_DIR=$(pwd)/Software
    mkdir -p $INSTALL_DIR/bin
    cd $INSTALL_DIR/bin
        
    #  Install nextflow (latest)
    wget -qO- https://get.nextflow.io | bash
    
    #  Bismark (0.23.0)
    
    wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz && \
        tar -xzf 0.23.0.tar.gz
        cp bismark* ../bin/
        cp coverage2cytosine ../bin/
        cd $INSTALL_DIR/
    
    #  fastqc (0.11.8)  -------------------------------------------------------
    
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
        unzip fastqc_v0.11.8.zip && \
        chmod -R 775 FastQC
        cp FastQC/fastqc bin/
        
    #  MethylDackel
    
    ##  Install libBigWig, a dependency
    git clone git@github.com:dpryan79/libBigWig.git && \
            cd libBigWig && \
            make prefix=$INSTALL_DIR install
            cd $INSTALL_DIR
      
    ##  Install htslib (1.10.2), a dependency     
    wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 -O htslib.tar.bz2 && \
            tar -xjf htslib.tar.bz2 && \
            cd htslib-1.10.2 && \
            ./configure prefix=$INSTALL_DIR && \
            make && \
            make install
            cd $INSTALL_DIR
            mv htslib-1.10.2 htslib
    
    ##  MethylDackel itself
    git clone https://github.com/dpryan79/MethylDackel.git && \
        cd MethylDackel && \
        make install CFLAGS="-O3 -Wall -I$INSTALL_DIR/include " LIBS="-L$INSTALL_DIR/lib" prefix=$INSTALL_DIR/bin LIBBIGWIG="$INSTALL_DIR/libBigWig/libBigWig.a"
        cd $INSTALL_DIR
        
    #  R (3.6.1) --------------------------------------------------------------
        
    #  Install R
    wget http://cran.rstudio.com/src/base/R-3/R-3.6.1.tar.gz && \
        tar -xf R-3.6.1.tar.gz && \
        cd R-3.6.1 && \
        ./configure --prefix=$INSTALL_DIR && \
        make && \
        make install
        cd $INSTALL_DIR
      
    #  Install packages that will be used by the pipeline
    ./R-3.6.1/bin/Rscript ../scripts/install_R_packages.R
    
    #  samblaster
    
    wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz && \
        tar -xz samblaster-v.0.1.26.tar.gz
    
    #  samtools (1.10)  -------------------------------------------------------
    
    wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2 && \
        tar -xjf samtools.tar.bz2 && \
        cd samtools-1.10 && \
        ./configure prefix=$INSTALL_DIR && \
        make && \
        make install
        cd $INSTALL_DIR
        
    #  Trim Galore!
        
    #  Fix any strict permissions which would not allow sharing software
    #  (and therefore the pipeline as a whole) with those in a group
    chmod 775 -R $INSTALL_DIR
    
else #  Java could not be found on the system
    echo "A java runtime could not be found or accessed. Is it installed and on the PATH? You can install it by running 'apt install default-jre', which requires sudo/ root privileges."
    echo "After installing Java, rerun this script to finish the installation procedure."
fi