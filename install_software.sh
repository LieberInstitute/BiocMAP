#!/bin/bash

#  IMPORTANT: run this script inside the base repository directory to ensure 
#             software locations are properly linked.
#
#  Usage:  bash install_software.sh [installation_type]
#
#    installation_type may be "docker", "local", "conda", "jhpce", or
#    "singularity":
#        docker:      user plans to run BiocMAP with docker to manage software
#                     dependencies
#        local:       user wishes to install all dependencies locally (regardless
#                     of whether BiocMAP will be run on a cluster or with
#                     local resources)
#        conda:       required software is installed within a conda environment
#        jhpce:       user is setting up BiocMAP on the JHPCE cluster
#        singularity: user plans to run BiocMAP with singularity to manage
#                     software dependencies

set -e

REPO_NAME="BiocMAP"

if [ "$(basename $(pwd))" == "$REPO_NAME" ]; then

    echo "[BiocMAP] Making sure the repository is clean and ready for installation..."
    
    #rm -r Software
    rm -f test/*/*/samples.manifest
    rm -f test/*/*/rules.txt
    git checkout run_*_half_*.sh
    git checkout nextflow.config
    git checkout conf/*_half_*.config
    
else

    echo "[BiocMAP] Please only invoke the script from directly inside the '$REPO_NAME' directory!"
    exit 1
    
fi

if [[ "$1" == "docker" || "$1" == "singularity" ]]; then

    #  This is the docker image to be used for execution of R via docker/
    #  singularity
    R_container="libddocker/bioc_kallisto:3.17"
    
    #  Point to original repo's main script to facilitate pipeline sharing
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_*_half_*.sh
    
    #  Add docker/singularity configuration to each config profile in
    #  'nextflow.config'
    sed -i "s|includeConfig 'conf/\(.*\)_half_\(.*\)\.config'|includeConfig 'conf/\1_half_\2.config'\n        includeConfig 'conf/\1_half_$1.config'|" nextflow.config
    
    #  Remove a variable from configs that disables docker/singularity settings
    sed -i "/using_containers = false/d" conf/second_half_{jhpce,local,sge,slurm}.config
    
    BASE_DIR=$(pwd)
    mkdir -p $BASE_DIR/Software/bin
    cd $BASE_DIR/Software/bin
    
    #  Install nextflow (latest)
    echo "[BiocMAP] Installing nextflow..."
    wget -qO- https://get.nextflow.io | bash
        
    ###########################################################################
    #  Create the 'samples.manifest' and 'rules.txt' files for test samples
    ###########################################################################
    
    if [[ "$1" == "docker" ]]; then
        if [[ "$2" == "sudo" ]]; then
            command="sudo docker run"
        else
            command="docker run"
        fi
        
        echo "[BiocMAP] Setting up test files..."

        $command \
            -it \
            -u $(id -u):$(id -g) \
            -v $BASE_DIR/scripts:/usr/local/src/scripts/ \
            -v $BASE_DIR/test:/usr/local/src/test \
            $R_container \
            Rscript /usr/local/src/scripts/prepare_test_files.R -d $BASE_DIR
    else # using singularity
        echo "[BiocMAP] Pulling docker images and converting to singularity images..."
        cd $BASE_DIR
        
        #  Pull images in advance, since it seems to use very large amounts of
        #  memory to build the '.sif' file from each docker image (we don't
        #  want to allocate large amounts of memory in each process just for
        #  this purpose)
        mkdir -p docker/singularity_cache
        images=$(grep 'container = ' conf/*_half_singularity.config | tr -d " |'" | cut -d '=' -f 2 | sort -u)
        for image in $images; do
            echo "[BiocMAP] Pulling and converting ${image}..."

            image_name=$(echo $image | sed 's/[:\/]/-/g').sif
            singularity pull docker/singularity_cache/$image_name docker://$image
        done

        echo "[BiocMAP] Setting up test files..."
        
        singularity exec \
            -B $BASE_DIR/scripts:/usr/local/src/scripts/ \
            -B $BASE_DIR/test:/usr/local/src/test \
            docker://$R_container \
            Rscript /usr/local/src/scripts/prepare_test_files.R -d $BASE_DIR
        
        #  Set modules used correctly for JHPCE users
        sed -i "/module = '.*\/.*'/d" conf/*_half_jhpce.config
        sed -i "s|cache = 'lenient'|cache = 'lenient'\n    module = 'singularity/3.6.0'|" conf/*_half_jhpce.config
        sed -i "s|module load nextflow|module load nextflow\nmodule load singularity/3.6.0|" run_*_half_jhpce.sh
    fi
        
    echo "[BiocMAP] Done."
    
elif [ "$1" == "jhpce" ]; then

    echo "[BiocMAP] User selected set-up at JHPCE. Installing any missing R packages..."
    module load conda_R/4.3
    Rscript scripts/install_r_packages.R
    
    echo "[BiocMAP] Setting up test files..."
    Rscript scripts/prepare_test_files.R -d $(pwd)
    
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_first_half_jhpce.sh
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_second_half_jhpce.sh
    
    echo "[BiocMAP] Done."
    
elif [ "$1" == "conda" ]; then
    
    BASE_DIR=$(pwd)
    mkdir -p $BASE_DIR/Software/bin
    cd $BASE_DIR/Software/bin
    
    #  Install nextflow (latest)
    echo "[BiocMAP] Installing nextflow..."
    wget -qO- https://get.nextflow.io | bash
    
    cd $BASE_DIR

    #  Create an initial small conda environment, just containing mamba
    echo "[BiocMAP] Creating a conda environment containing required software..."
    source $(conda info --base)/etc/profile.d/conda.sh
    conda env create -p $PWD/conda/pipeline_env -f conda/mamba.yml
    conda activate $PWD/conda/pipeline_env
    
    #  Install software using mamba rather than conda
    mamba install -y -c bioconda -c conda-forge r-essentials=4.3.1 r-base=4.3.1 bismark=0.23.0 fastqc=0.11.8 kallisto=0.46.1 methyldackel=0.6.0 samblaster=0.1.26 samtools=1.12 trim-galore=0.6.6
    
    #  Install Bioc R packages using mamba
    mamba install -y -c bioconda -c conda-forge bioconductor-bsseq=1.36.0 bioconductor-genomicranges=1.52.0 bioconductor-hdf5array=1.28.1 bioconductor-biocparallel=1.34.2
    
    #  Install remaining non-Bioc packages
    Rscript scripts/install_r_packages.R
    
    #   If an NVIDIA GPU is available (if 'nvidia-smi' works as a command), install Arioc
    if nvidia-smi ; then
        echo "[BiocMAP] Installing Arioc, which isn't available as a conda package..."

        cd $BASE_DIR/Software/
        mkdir arioc
        cd arioc
            
        wget https://github.com/RWilton/Arioc/releases/download/v1.43/Arioc.x.143.zip
        unzip Arioc.x.143.zip
        
        cd src
        make clean
        make AriocE
        make AriocU
        make AriocP
        
        cd $BASE_DIR
        cp Software/arioc/bin/* conda/pipeline_env/bin/
    else
        echo "[BiocMAP] Warning: skipping Arioc installation, as no NVIDIA GPU was detected."
    fi
    
    echo "[BiocMAP] Setting up test files..."
    Rscript scripts/prepare_test_files.R -d $(pwd)
    conda deactivate
    
    echo "[BiocMAP] Configuring main and config files..."
    
    #  Point to original repo's main script to facilitate pipeline sharing
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_*_half_*.sh
    
    #  Add 'conda' specification to generic process scope; remove use of
    #  modules in JHPCE configs
    sed -i "s|cache = 'lenient'|cache = 'lenient'\n    conda = '$PWD/conda/pipeline_env'|" conf/*_half_*.config
    sed -i "/module = '.*\/.*'/d" conf/*_half_jhpce.config
    
    echo "[BiocMAP] Done."
    
elif [ "$1" == "local" ]; then
    #  Verify java can be executed, since this is a pre-requisite
    if [ -x "$(command -v java)" ]; then
        echo "[BiocMAP] Found a java runtime. Proceeding with the setup..."
        
        #  Point to original repo's main script to facilitate pipeline sharing
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_*_half_*.sh
        
        #  Add some configuration of environment variables for local runs
        sed -i '5s|\(.*\)|\1\n//  Add locally installed software to path-related environment variables\nrepoDir=System.getProperty("user.dir")\nenv.LD_LIBRARY_PATH="\$repoDir/Software/lib:System.getenv('"'"'LD_LIBRARY_PATH'"'"')"\nenv.PATH="\$repoDir/Software/bin:System.getenv('"'"'PATH'"'"')"\n|' conf/*_half_local.config
        
        BASE_DIR=$(pwd)
        INSTALL_DIR=$(pwd)/Software
        mkdir -p $INSTALL_DIR/bin
        cd $INSTALL_DIR/bin
            
        #  Install nextflow (latest)
        wget -qO- https://get.nextflow.io | bash
        cd $INSTALL_DIR
        
        #  Arioc (1.43) -------------------------------------------------------
        
        #   If an NVIDIA GPU is available (if 'nvidia-smi' works as a command), install Arioc
        if nvidia-smi ; then
            mkdir arioc
            cd arioc
            
            wget https://github.com/RWilton/Arioc/releases/download/v1.43/Arioc.x.143.zip
            unzip Arioc.x.143.zip
            
            cd src
            make clean
            make AriocE
            make AriocU
            make AriocP
            
            cd $INSTALL_DIR
            cp arioc/bin/* bin/
        else
            echo "[BiocMAP] Warning: skipping Arioc installation, as no NVIDIA GPU was detected."
        fi
        
        #  Bismark (0.23.0) ---------------------------------------------------
        
        wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz
        tar -xzf 0.23.0.tar.gz
        cp Bismark-0.23.0/bismark* bin/
        cp Bismark-0.23.0/coverage2cytosine bin/
        
        #  cmake, which kallisto need to build --------------------------------
        git clone https://gitlab.kitware.com/cmake/cmake.git
        cd cmake
        ./bootstrap --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        
        #  fastqc (0.11.8)  ---------------------------------------------------
        
        wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
        unzip fastqc_v0.11.8.zip
        chmod -R 775 FastQC
        cp FastQC/fastqc bin/
        
        #  kallisto (0.46.1)  -------------------------------------------------
                    
        wget https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
        tar -xzf v0.46.1.tar.gz
        cd kallisto-0.46.1
        mkdir build
        cd build
        $INSTALL_DIR/bin/cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
        make
        make prefix=$INSTALL_DIR install
        cd $INSTALL_DIR
            
        #  MethylDackel (0.5.2) -----------------------------------------------
        
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
          
        #  Install packages and set up test files
        echo "[BiocMAP] Installing R packages..."
        Rscript ../scripts/install_r_packages.R
        
        echo "[BiocMAP] Setting up test files..."
        Rscript ../scripts/prepare_test_files.R -d $BASE_DIR
        
        #  samblaster (v.0.1.26) ----------------------------------------------
        
        wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz
        tar -xzf samblaster-v.0.1.26.tar.gz
        cd samblaster-v.0.1.26
        make
        cp samblaster ../bin/
        cd $INSTALL_DIR
        
        #  samtools (1.10)  ---------------------------------------------------
        
        wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2
        tar -xjf samtools.tar.bz2
        cd samtools-1.10
        ./configure prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
            
        #  Trim Galore! (0.6.6) -----------------------------------------------
        
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
        echo "[BiocMAP] A java runtime could not be found or accessed. Is it installed and on the PATH? You can install it by running 'apt install default-jre', which requires sudo/ root privileges."
        echo "[BiocMAP] After installing Java, rerun this script to finish the installation procedure."
    fi
else # neither "docker", "local", "conda", "jhpce", nor "singularity" were chosen
    
    echo '[BiocMAP] Error: please specify "docker", "local", "conda", "jhpce", or "singularity" and rerun this script.'
    echo '[BiocMAP]     eg. bash install_software.sh "local"'
    exit 1
    
fi
