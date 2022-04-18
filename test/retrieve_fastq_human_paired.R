library('synapser')

dest_dir = 'human/paired/fastq'
num_reads = 100000

syn_ids = c(
    'syn7091230', 'syn7092119', 'syn7092951', 'syn7093942', 'syn7096852',
    'syn7097424'
)

#   Enter credentials manually
synLogin(email = "", authToken = "")

#   Pull files individually
sapply(syn_ids, function(x) synGet(x, downloadLocation = 'temp'))

#   Subset to a small number of reads
fastq_temp = list.files('temp', full.names = TRUE)
for (fastq in fastq_temp) {
    command = paste0(
        'gunzip -c ', fastq, ' | head -n ', num_reads, ' | gzip > ',
        file.path(dest_dir, basename(fastq))
    )
    
    system(command)
}
