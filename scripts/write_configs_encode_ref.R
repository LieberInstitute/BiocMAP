library('getopt')

spec = matrix(
    c(
        'ref', 'r', 1, 'character', 'reference genome',
        'gap_seed', 'g', 1, 'character', 'gapped seed',
        'nongap_seed', 'n', 1, 'character', 'non-gapped seed'
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

seq_names = readLines(list.files(pattern='chr_names_.*'))

for (seed in c(opt$gap_seed, opt$nongap_seed)) {
    #  The beginning portion of the config. Here the literal
    #  '[future_work_dir]' is written and will be replaced with the working
    #  directory used in the EncodeReference process
    config_lines = c('<?xml version="1.0" encoding="utf-8"?>',
                     paste0('<AriocE seed="', seed, '" gpuMask="0x00000001" maxJ="200">'),
                     paste0('  <dataIn sequenceType="R" srcId="0" filePath="[future_work_dir]">'))
    
    #  The sequences to encode
    seq_lines = sapply(1:length(seq_names), function(i) {
                          paste0('    <file SN="', seq_names[i], '" subId="', i, '">', seq_names[i], '.fa</file>')
                      })
                      
    #  The remaining portion, again using the literal '[future_work_dir]' to be
    #  replaced later
    config_lines = c(config_lines, seq_lines,
                     '  </dataIn>',
                     '  <dataOut>',
                     paste0('    <path>[future_work_dir]</path>'),
                     '  </dataOut>',
                     '</AriocE>')
    
    filename = paste0('encode_ref_', seed, '.cfg')
    writeLines(config_lines, filename)
}
