library('getopt')

spec <- matrix(c('ref', 'r', 1, 'character', 'reference genome',
                 'dir', 'd', 1, 'character', 'output directory'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

seq_names = readLines(list.files(pattern='chr_names_.*'))

for (seed in c("hsi25_0_32_CT", "ssi84_2_30_CT")) {
    #  The beginning portion of the config
    config_lines = c('<?xml version="1.0" encoding="utf-8"?>',
                     paste0('<AriocE seed="', seed, '" maxJ="200">'),
                     paste0('  <dataIn sequenceType="R" srcId="0" filePath="', getwd(), '">'))
    
    #  The sequences to encode
    seq_lines = sapply(1:length(seq_names), function(i) {
                          paste0('    <file SN="', seq_names[i], '" subId="', i, '">', seq_names[i], '.fa</file>')
                      })
                      
    #  The remaining portion
    config_lines = c(config_lines, seq_lines,
                     '  </dataIn>',
                     '  <dataOut>',
                     paste0('    <path>', opt$dir, '</path>'),
                     '  </dataOut>',
                     '</AriocE>')
    
    #  Write config for this seed type
    if (seed == "hsi25_0_32_CT") {
        filename = 'encode_ref_gap.cfg'
    } else {
        filename = 'encode_ref_nongap.cfg'
    }

    writeLines(config_lines, file.path(getwd(), filename))
}
