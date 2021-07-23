#!/bin/bash
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G
#$ -o ./run_first_half_jhpce.log
#$ -e ./run_first_half_jhpce.log
#$ -cwd

#  After running 'install_software.sh', this should point to the directory
#  where SPEAQeasy was installed, and not say "$PWD"
ORIG_DIR=$PWD

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $ORIG_DIR/first_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "single" \
    --reference "hg38" \
    --trim_mode "force" \
    -profile first_half_jhpce
