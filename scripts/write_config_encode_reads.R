library('getopt')

spec <- matrix(c('paired', 'p', 1, 'character', '"single" or "paired"',
                 'writeDir', 'd', 1, 'character', 'directory to put encoded reads in',
                 'prefix', 'x', 1, 'character', 'prefix to uniquely identify sample'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

if (opt$paired == "paired") {
    id = strsplit(system('ls *.fastq', intern=TRUE)[1], "_1.fastq", fixed=TRUE)[[1]]
} else {
    id = strsplit(system('ls *.fastq', intern=TRUE), ".fastq", fixed=TRUE)[[1]]
}

#  The processed manifest outputted from the "Manifest" process will be in
#  the working directory. Use this to assign this particular sample an integer ID
man = read.table('arioc_samples.manifest', sep = ' ', header = FALSE, stringsAsFactors = FALSE)
idRowNum = match(id, man[,ncol(man)])
idNum = as.integer(idRowNum / 2) + 1

#############################################################
#    Gather lines into a vector for writing
#############################################################

#  Starting lines
config_lines = c('<?xml version="1.0" encoding="utf-8"?>', '',
                 '<AriocE>',
                 '  <dataIn sequenceType="Q">')
                 
#  Lines related to the particular reads used for this alignment
if (opt$paired == "paired") {
    config_lines = c(config_lines,
                     paste0('    <file subId="', idNum, '" mate="1">', man[idRowNum, 5], '_1.fastq</file>'),
                     paste0('    <file subId="', idNum, '" mate="2">', man[idRowNum, 5], '_2.fastq</file>'))
} else {
    config_lines = c(config_lines,
                     paste0('    <file subId="', idNum, '">', man[idRowNum, 3], '.fastq</file>'))
}
config_lines = c(config_lines, '  </dataIn>')

#  Output
config_lines = c(config_lines,
                 '  <dataOut>',
                 paste0('    <path>', opt$writeDir, '</path>'),
                 '  </dataOut>',
                 '</AriocE>')

#############################################################
#    Write the config to file
#############################################################

writeLines(config_lines, paste0(opt$prefix, '_encode_reads.cfg'))