#!/bin/bash
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=50G
#$ -o ./run_jhpce_qsub.log
#$ -e ./run_jhpce_qsub.log
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow main.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile jhpce \
    -resume
