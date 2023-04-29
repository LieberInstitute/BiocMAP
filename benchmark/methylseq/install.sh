module load nextflow/22.10.7

#   Prompt pull of pipeline
nextflow run nf-core/methylseq

#   Convert the 'samples.manifest' into a nf-core-compatible format
module load conda_R/4.2.x
Rscript -e "
a = read.table(
    '../samples.manifest',
    col.names = c('fastq_1', 'md1', 'fastq_2', 'md2', 'sample')
)
write.csv(
    a[, c('sample', 'fastq_1', 'fastq_2')], '../sample_sheet.csv', row.names = FALSE,
    quote = FALSE
)
"
