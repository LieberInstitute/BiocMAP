#!/bin/bash
#SBATCH -t 30:00
#SBATCH --mem=10G

nextflow=/scratch/groups/ajaffe1/rna_sp/RNAsp/Software/nextflow
module load java

$nextflow run main.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile marcc
