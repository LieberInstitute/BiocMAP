#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=6
#SBATCH -N 1

cd /scratch/groups/ajaffe1/WGBS-Pipeline
nextflow=/scratch/groups/ajaffe1/rna_sp/RNAsp/Software/nextflow
module load java

$nextflow run main.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile marcc \
    -resume
