#!/bin/bash
#$ -cwd
#$ -o render.log
#$ -e render.log
#$ -l bluejay,mf=30G,h_vmem=30G

RMD_FILE=$1

module load conda_R/4.1.x
module load tex/2017

#  Style code
Rscript -e "styler::style_file('$RMD_FILE', transformers = biocthis::bioc_style())"

#  Render into PDF (as defined in the YAML header)
Rscript -e "rmarkdown::render('$RMD_FILE', clean = FALSE)"
