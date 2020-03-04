library('getopt')

spec <- matrix(c('paired', 'p', 1, 'character', '"single" or "paired"',
                 'writeDir', 'd', 1, 'character', 'directory to put encoded reads in',
                 'prefix', 'x', 1, 'character', 'prefix to uniquely identify sample'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

print("Writing AriocE configs...")

if (opt$paired == "paired") {
    id = strsplit(system('ls *.fastq', intern=TRUE)[1], "_1.fastq", fixed=TRUE)[[1]]
    
    ###################################################################################
    #  Downstream utilities like Bismark and samblaster are incapable of recognizing
    #  Illumina Casava-style FASTQ read-IDs (cannot distinguish mate ID from read-ID).
    #  This code identifies the format style of the FASTQ read ids for a pair of
    #  mates, and forwards information about QNAME and read group to use for a given
    #  FASTQ sequence to AriocE.
    ###################################################################################
    
    id1 = system(paste0('cat ', id, '_1.fastq | grep "@" | head -n 1'), intern=TRUE)
    id2 = system(paste0('cat ', id, '_2.fastq | grep "@" | head -n 1'), intern=TRUE)
    
    #----------------------------------------------------------------------------------
    #  Tries to recognize the read-ID style (checks "old illumina", "illumina casava")
    #----------------------------------------------------------------------------------
    
    print("Examining read ID format in the FASTQ files since data is paired...")
    recognized = FALSE
    
    temp = strsplit(id1, ':')[[1]]
    cond1 = length(temp) == 5
    cond2 = all(c('#', '/') %in% strsplit(temp[length(temp)], ''))
    if (cond1 & cond2) {
        recognized = TRUE
        
        #  Make sure the read IDs for the first site are the same between reads
        #  (aside from mate ID)
        first_id = temp[1:4]
        second_id = strsplit(id2, ':')[[1]][1:4]
        if (!all(first_id == second_id)) {
            stop("Read IDs for the first site differ by more than mate ID, between reads")
        }
        
        #  Make sure mate ID differs
        first_id = temp[5]
        second_id = strsplit(id2, ':')[[1]][5]
        if (first_id == second_id) {
            stop("Despite old Illumina-style read-ID, all fields (including mate ID) were equivalent between reads.")
        }
        
        print("FASTQ read-id is in the old Illumina format, and looks as expected.")
        
        qname_field = ''
        
        #  Use flowcell id as the "read group"
        rg_line = '    <rg ID="*:(*):" PL="ILLUMINA" />'
    }
    
    cond1 = length(temp) == 10
    cond2 = ' ' %in% strsplit(temp[7], '')[[1]]
    if (cond1 & cond2) {
        recognized = TRUE
        
        #  Make sure the read IDs for the first site are the same between reads
        #  (aside from mate ID)
        first_id = temp[c(1:6, 8:10)]
        second_id = strsplit(id2, ':')[[1]][c(1:6, 8:10)]
        if (!all(first_id == second_id)) {
            stop("Read IDs for the first site differ by more than mate ID, between reads")
        }
        
        #  Make sure mate ID differs
        first_id = temp[7]
        second_id = strsplit(id2, ':')[[1]][7]
        if (first_id == second_id) {
            stop("Despite Illumina casava-style read-ID, all fields (including mate ID) were equivalent between reads.")
        }
        
        print("FASTQ read-id is in the Illumina casava format, and looks as expected.")
        
        #  Many utilites fail to recognize read IDs as the same in this format (don't know where mate ID is),
        #  so we signal here to AriocE to remove the mate ID from the output SAM's QNAME field
        qname_field = 'QNAME="*:*:*:(*:*:*:*) "'
        
        #  Use flowcell id and lane as the "read group"
        rg_line = '    <rg ID="*:*:(*:*):" PL="ILLUMINA" />'
    }
    
    if (!recognized) {
        mess = "Warning: could not recognize the FASTQ read-id format. Checked Illumina (old and casava-style)."
        mess = paste(mess, "Downstream utilities like samblaster and BME might not recognize mate identifier, which")
        mess = paste(mess, "would halt analysis.")
        
        print(mess)
    }
    
} else {
    id = strsplit(system('ls *.fastq', intern=TRUE), ".fastq", fixed=TRUE)[[1]]
    qname_field = ''
    rg_line = ''
}

#  The processed manifest outputted from the "Manifest" process will be in
#  the working directory. Use this to assign this particular sample an integer ID
man = read.table('arioc_samples.manifest', sep = ' ', header = FALSE, stringsAsFactors = FALSE)
idRowNum = match(opt$prefix, man[,ncol(man)])
idNum = as.integer(idRowNum / 2) + 1

#############################################################
#    Gather lines into a vector for writing
#############################################################

#  Starting lines
config_lines = c('<?xml version="1.0" encoding="utf-8"?>', '',
                 '<AriocE>',
                 paste0('  <dataIn sequenceType="Q" ', qname_field, '>'),
                 rg_line)
                 
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

print("Done writing AriocE configs.")
