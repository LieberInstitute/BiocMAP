#!/bin/bash
#SBATCH --output=run_second_half_slurm.log
#SBATCH --mem=25G

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/bin/nextflow second_half.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile second_half_slurm
