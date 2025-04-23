#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=BiocMAP_second_module
#SBATCH -o ./run_second_half_jhpce.log
#SBATCH -e ./run_second_half_jhpce.log

#  After running 'install_software.sh', this should point to the directory
#  where BiocMAP was installed, and not say "$PWD"
ORIG_DIR=$PWD

module load nextflow/20.01.0
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $ORIG_DIR/second_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "paired" \
    --reference "hg38" \
    -profile second_half_jhpce

#   Log successful runs on non-test data in a central location. Please adjust
#   the log path here if it is changed at the top!
bash $ORIG_DIR/scripts/track_runs.sh $PWD/run_second_half_jhpce.log second
