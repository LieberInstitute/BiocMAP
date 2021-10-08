library('jaffelab')
library('getopt')

spec = matrix(
    c("repo_dir", "d", 1, "character", "path to this repo"),
    byrow = TRUE,
    ncol = 5
)
opt = getopt(spec)


for (pairing in c('single', 'paired')) {
    for (species in c('human', 'mouse')) {
        base_dir = file.path('test', species, pairing)
        
        if (pairing == 'paired') {
            #  Create samples.manifest
            r1 = file.path(opt$repo_dir, list.files(file.path(base_dir, 'fastq'), pattern = '.*_1\\.fastq\\.gz$', full.names = TRUE))
            r2 = file.path(opt$repo_dir, list.files(file.path(base_dir, 'fastq'), pattern = '.*_2\\.fastq\\.gz$', full.names = TRUE))
            ids = ss(basename(r1), '_1')
            
            man = paste(r1, 0, r2, 0, ids, sep='\t')
            writeLines(man, con = file.path(base_dir, 'samples.manifest'))
            
            #  Create rules.txt
            rules = c(
                paste0('manifest = ', file.path(opt$repo_dir, base_dir, 'samples.manifest')),
                paste0('sam = ', file.path(opt$repo_dir, base_dir, 'bam', '[id].cfus.bam*')),
                paste0('arioc_log = ', file.path(opt$repo_dir, base_dir, 'logs', '[id]_alignment.log')),
                paste0('trim_report = ', file.path(opt$repo_dir, base_dir, 'logs', '[id]_was_trimmed.log')),
                paste0('fastqc_log_last = ', file.path(opt$repo_dir, base_dir, 'logs', '[id]_[12]_fastqc', 'summary.txt'))
            )
            writeLines(rules, con = file.path(base_dir, 'rules.txt'))
        } else {
            #  Create samples.manifest
            r1 = file.path(opt$repo_dir, list.files(file.path(base_dir, 'fastq'), pattern = '.*\\.fastq\\.gz$', full.names = TRUE))
            ids = ss(basename(r1), '\\.')
            
            man = paste(r1, 0, ids, sep='\t')
            writeLines(man, con = file.path(base_dir, 'samples.manifest'))
            
            #  Create rules.txt
            rules = c(
                paste0('manifest = ', file.path(opt$repo_dir, base_dir, 'samples.manifest')),
                paste0('sam = ', file.path(opt$repo_dir, base_dir, 'bam', '[id].mfus.bam*')),
                paste0('arioc_log = ', file.path(opt$repo_dir, base_dir, 'logs', '[id]_alignment.log')),
                paste0('trim_report = ', file.path(opt$repo_dir, base_dir, 'logs', '[id]_was_trimmed.log')),
                paste0('fastqc_log_last = ', file.path(opt$repo_dir, base_dir, 'logs', '[id]_fastqc', 'summary.txt'))
            )
            writeLines(rules, con = file.path(base_dir, 'rules.txt'))
        }
    }
}
