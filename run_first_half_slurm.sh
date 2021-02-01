#!/bin/bash
#SBATCH --output=run_first_half_slurm.log
#SBATCH --mem=25G

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/bin/nextflow first_half.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile first_half_slurm
