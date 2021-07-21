#!/bin/bash
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G
#$ -o ./run_first_half_jhpce.log
#$ -e ./run_first_half_jhpce.log
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow first_half.nf \
    --sample "paired" \
    --reference "hg38" \
    --trim_mode "force" \
    -profile first_half_jhpce
