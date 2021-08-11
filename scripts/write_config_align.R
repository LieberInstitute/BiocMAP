library('getopt')

spec <- matrix(c('paired', 'p', 1, 'character', '"single" or "paired"',
                 'prefix', 'f', 1, 'character', 'prefix to uniquely identify sample',
                 'allAlign', 'a', 1, 'logical', 'whether to output nonconcordant alignments',
                 'arioc_opts', 'o', 1, 'character', 'main alignment options',
                 'gapped_opts', 'g', 1, 'character', 'gapped seed options',
                 'nongapped_opts', 'n', 1, 'character', 'nongapped seed options',
                 'x_opts', 'x', 1, 'character', 'opts to pass to <X>',
                 'q_opts', 'q', 1, 'character', 'file-related options'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

print("Writing configs for AriocP/AriocU...")

if (opt$paired == "paired") {
  exec_name = 'AriocP'
  id = strsplit(list.files(pattern=".*_val_1\\.fq"), "_val_1.fq", fixed=TRUE)[[1]]
} else {
  exec_name = 'AriocU'
  id = strsplit(list.files(pattern='.*_trimmed\\.fq'), "_trimmed.fq", fixed=TRUE)[[1]]
}

#  The processed manifest outputted from the "Merging" process will be in
#  the working directory. Infer the directory containing the input FASTQs from
#  the first sample in the manifest (thus the assumption is made: the manifest
#  paths all match the input directory, params.input)
man = read.table('arioc_samples.manifest', sep = ' ', header = FALSE, stringsAsFactors = FALSE)
idRowNum = match(opt$prefix, man[,ncol(man)])
idNum = as.integer(idRowNum / 2) + 1

if (opt$allAlign) {
    if (opt$paired == "paired") {
        sam_outputs = c("c", "r", "d", "u")
    } else {
        sam_outputs = c("m", "u")
    }
} else {
    if (opt$paired == "paired") {
        sam_outputs = "c"
    } else {
        sam_outputs = "m"
    }
}

#############################################################
#    Gather lines into a vector for writing
#############################################################

#  Lines related to alignment settings: here the literal '[future_work_dir]' is
#  written and will be replaced with the working directory used in the
#  AlignReads process
config_lines = c(
    '<?xml version="1.0" encoding="utf-8"?>',
    opt$arioc_opts,
    '  <R>[future_work_dir]</R>',
    '',
    opt$nongapped_opts,
    opt$gapped_opts,
    opt$x_opts, 
    '',
    '  <Q filePath="[future_work_dir]/encoded_reads">'
)

#  Lines related to the particular reads used for this alignment
if (opt$paired == "paired") {
    config_lines = c(config_lines,
                     paste0('    <paired subId="', idNum, '">'),
                     paste0('      <file>', man[idRowNum, 5], '_val_1</file>'),
                     paste0('      <file>', man[idRowNum, 5], '_val_2</file>'),
                     '    </paired>',
                     '  </Q>',
                     '',
                     '  <A overwrite="true" pairOrientation="c" pairCollision="ocd" pairFragmentLength="500" cigarFormat="MID" mapqVersion="2">')
} else {    
    config_lines = c(config_lines,
                     paste0('    <unpaired subId="', idNum, '">'),
                     paste0('      <file>', man[idRowNum, 3], '_trimmed</file>'),
                     '    </unpaired>',
                     '  </Q>',
                     '',
                     '  <A overwrite="true" cigarFormat="MID" mapqVersion="2">')
}

#  Lines related to sam output format(s)
for (out_type in sam_outputs) {
    config_lines = c(config_lines,
                     paste0('    <sam report="', out_type, '">./</sam>'))
}

config_lines = c(config_lines, 
                 '  </A>',
                 paste0('</', exec_name, '>'))

#############################################################
#    Write the config to file
#############################################################

writeLines(config_lines, paste0(opt$prefix, '_align_reads.cfg'))

print("Done writing all configs for Arioc.")
