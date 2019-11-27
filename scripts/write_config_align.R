library('getopt')

spec <- matrix(c('paired', 'p', 1, 'character', '"single" or "paired"',
                 'prefix', 'x', 1, 'character', 'prefix to uniquely identify sample',
                 'batchSize', 'b', 1, 'character', 'GPU batch size for Arioc',
                 'repoDir', 'd', 1, 'character', 'base path to repo',
                 'ref', 'r', 1, 'character', 'reference genome',
                 'allAlign', 'a', 1, 'logical', 'whether to output nonconcordant alignments'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

if (opt$paired == "paired") {
  exec_name = 'AriocP'
  id = strsplit(system('ls *.fastq', intern=TRUE)[1], "_1.fastq", fixed=TRUE)[[1]]
} else {
  exec_name = 'AriocU'
  id = strsplit(system('ls *.fastq', intern=TRUE), ".fastq", fixed=TRUE)[[1]]
}

#  The processed manifest outputted from the "Manifest" process will be in
#  the working directory. Infer the directory containing the input FASTQs from
#  the first sample in the manifest (thus the assumption is made: the manifest
#  paths all match the input directory, params.input)
man = read.table('arioc_samples.manifest', sep = ' ', header = FALSE, stringsAsFactors = FALSE)
idRowNum = match(id, man[,ncol(man)])
idNum = as.integer(idRowNum / 2) + 1

if (opt$allAlign) {
    sam_outputs = c("c", "r", "d", "u")
} else {
    sam_outputs = c("c")
}

#############################################################
#    Gather lines into a vector for writing
#############################################################

#  Lines related to alignment settings
config_lines = c('<?xml version="1.0" encoding="utf-8"?>',
                 paste0('<', exec_name, ' gpuMask="0x00000001" batchSize="', opt$batchSize, '" verboseMask="0xE0000007">'),
                 paste0('  <R>', opt$repoDir, '/ref/', opt$ref, '/encoded_ref</R>'), '',
                 '  <nongapped seed="ssi84_2_30_CT" maxJ="200" maxMismatches="5"/>',
                 '  <gapped seed="hsi25_0_32_CT" Wmxgs="2,6,5,3" Vt="L,0,1" maxJ="20" seedDepth="4"/>',
                 '  <X watchdogInterval="60" cgaReserved="24M" useHinGmem="1" useJinGmem="0" useHJinGPmem="0" serialLUTinit="1" />', '',
                 paste0('  <Q filePath="', opt$repoDir, '/Arioc/temp_encoded_reads">'))

#  Lines related to the particular reads used for this alignment
if (opt$paired == "paired") {
    config_lines = c(config_lines,
                     paste0('    <paired subId="', idNum, '">'),
                     paste0('      <file>', man[idRowNum, 5], '_R1</file>'),
                     paste0('      <file>', man[idRowNum, 5], '_R2</file>'),
                     '    </paired>')
} else {    
    config_lines = c(config_lines,
                     paste0('    <unpaired subId="', idNum, '">'),
                     paste0('      <file>', man[idRowNum, 3], '</file>'),
                     '    </unpaired>')
}

#  Lines related to sam output format(s)
config_lines = c(config_lines,
                 '  </Q>', '',
                 '  <A overwrite="true" pairOrientation="c" pairCollision="ocd" pairFragmentLength="500" cigarFormat="MID" mapqVersion="2">')
for (out_type in sam_outputs) {
    config_lines = c(config_lines,
                     paste0('    <sam report="', out_type, '">./</sam>'))
}

config_lines = c(config_lines, '  </A>', paste0('</', exec_name, '>'))

#############################################################
#    Write the config to file
#############################################################

writeLines(config_lines, paste0(opt$prefix, '_align_reads.cfg'))