#!/bin/bash
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G
#$ -o logs/run_second_half_jhpce_2.log
#$ -e logs/run_second_half_jhpce_2.log
#$ -cwd

REPO_DIR=$(git rev-parse --show-toplevel)
ANN_DIR=/dcs04/lieber/lcolladotor/BiocMAP_benchmark_LIBD001/BiocMAP_ref_temp
BASE_WORK_DIR=/dcs04/lieber/lcolladotor/BiocMAP_benchmark_LIBD001/BiocMAP

mkdir -p $BASE_WORK_DIR

module load nextflow/20.01.0
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $REPO_DIR/benchmark/BiocMAP_dir/BiocMAP/second_half.nf \
    --annotation "$ANN_DIR" \
    --sample "paired" \
    --reference "hg38" \
    --input "$BASE_WORK_DIR/out" \
    -w "$BASE_WORK_DIR/work" \
    --output "$BASE_WORK_DIR/out" \
    -with-report "logs/second_half_exec_report_2.html" \
    -profile second_half_jhpce
