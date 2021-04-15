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
    
    #  Bismark (0.23.0) -------------------------------------------------------
    
    ##  hisat2 (2.2.1), required for bismark_genome_preparation
    
    wget https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.tar.gz
    tar -xzf v2.2.1.tar.gz
    cd hisat2-2.2.1
    cp hisat2 $INSTALL_DIR/bin/
    cp hisat2-align* $INSTALL_DIR/bin/
    cp hisat2-build* $INSTALL_DIR/bin/
    cp hisat2-inspect* $INSTALL_DIR/bin/
    cp *.py $INSTALL_DIR/bin/
    cd $INSTALL_DIR
    
    wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz
    tar -xzf 0.23.0.tar.gz
    cp Bismark-0.23.0/bismark* bin/
    cp Bismark-0.23.0/coverage2cytosine bin/
    
    #  cmake, which kallisto need to build ------------------------------------
    git clone https://gitlab.kitware.com/cmake/cmake.git
    cd cmake
    ./bootstrap --prefix=$INSTALL_DIR
    make
    make install
    cd $INSTALL_DIR
    
    #  fastqc (0.11.8)  -------------------------------------------------------
    
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
    unzip fastqc_v0.11.8.zip
    chmod -R 775 FastQC
    cp FastQC/fastqc bin/
    
    #  kallisto (0.46.1)  -------------------------------------------------------------
                
    wget https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
    tar -xzf v0.46.1.tar.gz
    cd kallisto-0.46.1
    mkdir build
    cd build
    $INSTALL_DIR/bin/cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
    make
    make prefix=$INSTALL_DIR install
    cd $INSTALL_DIR
        
    #  MethylDackel (latest)
    
    ##  Install libBigWig (0.4.6), a dependency
    wget https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.6.tar.gz
    tar -xzf 0.4.6.tar.gz
    cd libBigWig-0.4.6
    make prefix=$INSTALL_DIR install
    cd $INSTALL_DIR
      
    ##  Install htslib (1.12), a dependency
    wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
    tar -xjf htslib-1.12.tar.bz2
    cd htslib-1.12
    ./configure
    make prefix=$INSTALL_DIR
    make prefix=$INSTALL_DIR install
    cd $INSTALL_DIR
    
    ##  MethylDackel itself (0.5.2)
    wget https://github.com/dpryan79/MethylDackel/archive/refs/tags/0.5.2.tar.gz
    tar -xzf 0.5.2.tar.gz
    cd MethylDackel-0.5.2
    make CFLAGS="-O3 -Wall -I${INSTALL_DIR}/include " LIBS="-L${INSTALL_DIR}/lib" LIBBIGWIG="${INSTALL_DIR}/lib/libBigWig.a"
    make install prefix=${INSTALL_DIR} CFLAGS="-O3 -Wall -I${INSTALL_DIR}/include " LIBS="-L${INSTALL_DIR}/lib" LIBBIGWIG="${INSTALL_DIR}/lib/libBigWig.a"
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
