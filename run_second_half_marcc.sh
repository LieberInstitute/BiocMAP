#!/bin/bash
#SBATCH -t 30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task 1
#SBATCH -N 1
#SBATCH --output=run_second_half_marcc.log

#  After running 'install_software.sh', this should point to the directory
#  where this repo was cloned, and not say "$PWD"
ORIG_DIR=$PWD

module load java
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

$ORIG_DIR/Software/bin/nextflow $ORIG_DIR/second_half.nf \
    --annotation "$HOME/scratch/nextflow_ref" \
    -w "$HOME/scratch/nextflow_work_second" \
    --output "$HOME/scratch/nextflow_out" \
    --sample "paired" \
    --reference "hg38" \
    -profile second_half_marcc
