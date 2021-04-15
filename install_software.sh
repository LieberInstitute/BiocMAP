#!/bin/bash

#  IMPORTANT: run this script inside the base repository directory to ensure 
#             software locations are properly linked.
#
#  Usage:  bash install_software.sh

set -e

#  Verify java can be executed, since this is a pre-requisite
if [ -x "$(command -v java)" ]; then
    echo "Found a java runtime. Proceeding with the setup..."
    
    INSTALL_DIR=$(pwd)/Software
    mkdir -p $INSTALL_DIR/bin
    cd $INSTALL_DIR/bin
        
    #  Install nextflow (latest)
    wget -qO- https://get.nextflow.io | bash
    cd $INSTALL_DIR
    
    #  Bismark (0.23.0)
    
    wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz
    tar -xzf 0.23.0.tar.gz
    cp Bismark-0.23.0/bismark* bin/
    cp Bismark-0.23.0/coverage2cytosine bin/
    
    #  fastqc (0.11.8)  -------------------------------------------------------
    
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
    unzip fastqc_v0.11.8.zip
    chmod -R 775 FastQC
    cp FastQC/fastqc bin/
        
    #  MethylDackel (latest)
    
    ##  Install libBigWig, a dependency
    git clone git@github.com:dpryan79/libBigWig.git
    cd libBigWig
    make prefix=$INSTALL_DIR install
    cd $INSTALL_DIR
      
    ##  Install htslib (1.10.2), a dependency     
    wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 -O htslib.tar.bz2
    tar -xjf htslib.tar.bz2
    cd htslib-1.10.2
    ./configure prefix=$INSTALL_DIR
    make
    make install
    cd $INSTALL_DIR
    mv htslib-1.10.2 htslib
    
    ##  MethylDackel itself (latest)
    git clone https://github.com/dpryan79/MethylDackel.git
    cd MethylDackel
    make install CFLAGS="-O3 -Wall -I$INSTALL_DIR/include " LIBS="-L$INSTALL_DIR/lib" prefix=$INSTALL_DIR/bin LIBBIGWIG="$INSTALL_DIR/libBigWig/libBigWig.a"
    cd $INSTALL_DIR
      
    #  Install packages that will be used by the pipeline
    Rscript ../scripts/install_R_packages.R
    
    #  samblaster (v.0.1.26)
    
    wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
    tar -xzf samblaster-v.0.1.26.tar.gz
    cd samblaster-v.0.1.26
    make
    cp samblaster ../bin/
    cd $INSTALL_DIR
    
    #  samtools (1.10)  -------------------------------------------------------
    
    wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2
    tar -xjf samtools.tar.bz2
    cd samtools-1.10
    ./configure prefix=$INSTALL_DIR
    make
    make install
    cd $INSTALL_DIR
        
    #  Trim Galore!
    
    wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz
    tar -xzf 0.6.6.tar.gz
    cp TrimGalore-0.6.6/trim_galore bin/
        
    #  Clean up compressed files
    rm $INSTALL_DIR/*.tar.gz
    rm $INSTALL_DIR/*.bz2
    rm $INSTALL_DIR/*.zip
        
    #  Fix any strict permissions which would not allow sharing software
    #  (and therefore the pipeline as a whole) with those in a group
    chmod 775 -R $INSTALL_DIR
    
else #  Java could not be found on the system
    echo "A java runtime could not be found or accessed. Is it installed and on the PATH? You can install it by running 'apt install default-jre', which requires sudo/ root privileges."
    echo "After installing Java, rerun this script to finish the installation procedure."
fi
