library('jaffelab')

man_dir = paste0(getwd(), '/test')
R1 = paste0(man_dir, '/', list.files(man_dir, pattern = '*_R1\\.fastq'))
R2 = paste0(man_dir, '/', list.files(man_dir, pattern = '*_R2\\.fastq'))
man_text = paste(R1, 0, R2, 0, basename(ss(R1, '_R1\\.fastq')), sep = "\t")

writeLines(man_text, paste0(man_dir, '/samples.manifest'))
