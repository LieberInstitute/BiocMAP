#!/bin/bash
#SBATCH -t 15:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task 3
#SBATCH -N 1
#SBATCH --output=run_first_half_marcc.log

#  After running 'install_software.sh', this should point to the directory
#  where this repo was cloned, and not say "$PWD"
ORIG_DIR=$PWD

module load java
export _JAVA_OPTIONS="-Xms3g -Xmx4g"

$ORIG_DIR/Software/bin/nextflow $ORIG_DIR/first_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "paired" \
    --reference "hg38" \
    -profile first_half_marcc
