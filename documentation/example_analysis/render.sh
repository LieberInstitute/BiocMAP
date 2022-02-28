#!/bin/bash
#$ -cwd
#$ -o render.log
#$ -e render.log
#$ -l mf=50G,h_vmem=50G

RMD_FILE=example_analysis.Rmd

module load conda_R/4.1.x
module load tex/2017

#  Style code
Rscript -e "styler::style_file('$RMD_FILE', transformers = biocthis::bioc_style())"

#  Render into PDF (as defined in the YAML header)
Rscript -e "rmarkdown::render('$RMD_FILE', clean = FALSE)"
