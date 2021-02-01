#!/bin/bash
#$ -l mem_free=25G,h_vmem=25G,h_fsize=800G
#$ -o ./run_second_half_sge.log
#$ -e ./run_second_half_sge.log
#$ -cwd

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/bin/nextflow second_half.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile second_half_sge
