#!/bin/bash
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G
#$ -o ./run_jhpce_qsub.log
#$ -e ./run_jhpce_qsub.log
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow main.nf \
    --sample "paired" \
    --reference "hg38" \
    --input "/dcl01/lieber/ajaffe/Nick/misc/manifests/WGBS_test_small" \
    --output "/dcl01/lieber/ajaffe/Nick/misc/manifests/WGBS_test_small/out" \
    -profile jhpce
