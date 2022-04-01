#!/bin/bash
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=800G
#$ -o ./run_second_half_jhpce.log
#$ -e ./run_second_half_jhpce.log
#$ -cwd

#  After running 'install_software.sh', this should point to the directory
#  where BiocMAP was installed, and not say "$PWD"
ORIG_DIR=$PWD

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $ORIG_DIR/second_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "paired" \
    --reference "hg38" \
    -profile second_half_jhpce

#   Log successful runs on non-test data in a central location. Please adjust
#   the log path here if it is changed at the top!
bash $ORIG_DIR/scripts/track_runs.sh $PWD/run_second_half_jhpce.log second
