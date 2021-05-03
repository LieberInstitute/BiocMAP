library('jaffelab')

get_value = function(rules, key, required) {
    this_line = rules[grep(key, ss(rules, '='))]
    
    return(gsub(' ', '', ss(this_line, '=', 2)))
}

#  Get sample IDs
rules = readLines('rules.txt')
rules = rules[-grep('#', rules)]

manifest = read.table(get_value(rules, 'manifest'), header = FALSE,
                      stringsAsFactors = FALSE)
ids = manifest[,ncol(manifest)]

######################################################
#  Arioc
######################################################

print("Extracting Arioc metrics...")

#  String patterns to remove via gsub for every value to put in a data frame:
#  note that patterns are removed in order of occurrence in this list
blacklist = c("(\\(.*\\))", "SAMrecords")
cleanseVec <- function(vec) {
    #  Remove blacklisted characters
    for (b in blacklist) {
        vec <- sapply(vec, function(x) gsub(b,"",x, perl=FALSE))
    }
    return(unname(vec))
}

col_names = c("pairs", "conc_pairs_total", "conc_pairs_1_mapping", "conc_pairs_many_mappings",
              "disc_pairs", "rejected_pairs", "unmapped_pairs","mates_not_in_paired_maps_total",
              "mates_NIPM_with_no_maps", "mates_NIPM_with_1_map", "mates_NIPM_with_many_maps",
              "total_mapped_mates", "duplicate_maps", "maxQlen", "max_diag_band_width", "TLEN mean")
              
f = paste0(ids, '_arioc.log')
stopifnot(file.exists(f[1]))

row_data = lapply(f, function(filename) {
    log_text = gsub(' ', '', system(paste0('cat ', filename, ' | cut -d " " -f 3-'), intern=TRUE))
    first_row = match("SAMoutput:", log_text) + 1
    log_text = log_text[first_row: (first_row + length(col_names) - 1)]
    
    as.numeric(cleanseVec(ss(log_text,":", 2)))
})

metrics = as.data.frame(matrix(unlist(row_data), nrow=length(row_data), byrow=TRUE))
colnames(metrics) = col_names

######################################################
#  BME logs
######################################################

filenames = paste0(ids, '_bme.log')
if (file.exists(filenames[1])) {
    print('Extracting BME metrics...')
    
    percs = lapply(filenames, function(f) {
        raw_chars = system(paste0('cat ', f, ' | grep "C methylated in C.. context:" | cut -d ":" -f 2'), intern=TRUE)
        return(as.numeric(sub("[:space:]|%", "", raw_chars)))
    })
    mat = matrix(unlist(percs), dimnames = list(ids, c("perc_M_CpG", "perc_M_CHG", "perc_M_CHH")), byrow=TRUE, ncol=3)
    
    metrics = cbind(metrics, mat)
} else {
    print('Skipping BME metrics (no files specified)...')
}

######################################################
#  Trimming Reports
######################################################
    
filepaths = paste0(ids, '_trim_report.log')

if (length(which(file.exists(filepaths))) > 0) {
    print('Extracting trimming metrics...')
    
    df_entries = c()
    
    is_first_iter = TRUE
    for (filename in filepaths) {
        if (file.exists(filename)) {
            a = readLines(filename)
            starts = which(a == '=== Summary ===')
            
            for (start_index in starts) {
                #  Subset to nonempty lines in the relevant section of the log
                b = a[(start_index+1):(start_index+12)]
                b = b[b != '' & b != "=== Adapter 1 ==="]
                
                #  Form appropriate column names from the lines
                if (is_first_iter) {
                    col_names = gsub('[ -]', '_', ss(b, ':'))
                    col_names = gsub('[()]', '', col_names)
                    is_first_iter = FALSE
                }
    
                #  Cut and format lines to extract values
                b_values = ss(b, ':', 2)
                b_values = ss(b_values, '[b;\\(]')
                b_values = gsub('[ ,]', '', b_values)
                
                #  Append to ongoing vector of values
                df_entries = c(df_entries, b_values)
            }
            
            #  The last column indicates that this sample was trimmed
            df_entries = c(df_entries, 'TRUE')
        } else {
            #  Fill in NAs for samples that weren't trimmed
            if (ncol(manifest) > 3) {
                col_len = 14
            } else {
                col_len = 7
            }
            
            #  The last column indicates that this sample wasn't trimmed
            df_entries = c(df_entries, rep(NA, col_len), 'FALSE')
        }
    }
    
    #  Use the correct colnames for paired-end reads
    if (ncol(manifest) > 3) {
        col_names = as.vector(outer(col_names, c('_R1', '_R2'), paste0))
    }
    
    #  Account for the last column indicating if the sample was trimmed
    col_names = c(col_names, 'was_trimmed')
    
    #  Form a data frame
    temp_df = matrix(df_entries, ncol=length(col_names), byrow=TRUE)
    temp_df = as.data.frame(temp_df)
    colnames(temp_df) = col_names
    
    #  Make the integer-containing columns numeric
    numeric_cols = which(!grepl('Sequence|was_trimmed', col_names))
    for (i in numeric_cols) {
        temp_df[,i] = as.numeric(temp_df[,i])
    }
    temp_df[,'was_trimmed'] = as.logical(temp_df[,'was_trimmed'])
    
    metrics = cbind(metrics, temp_df)
} else {
    print("Skipping trimming metrics (no samples appeared to be trimmed)...")
    metrics = cbind(metrics, data.frame('was_trimmed'=rep(FALSE, length(ids))))
}

######################################################
#  Lambda pseudoalignment
######################################################

f = paste0(ids, '_lambda_pseudo.log')
if (file.exists(f[1])) {
    print("Adding inferred BS-conversion efficiency from lambda alignment...")
    
    conv_eff = rep('empty', length(ids))
    for (i in 1:length(f)) {
        conv_eff[i] = system(paste0('tail -n 1 ', f[i], ' | cut -d ":" -f 2 | cut -d "%" -f 1 | tr -d " "'), intern=TRUE)
    }
    temp = colnames(metrics)
    metrics = cbind(metrics, as.numeric(conv_eff))
    colnames(metrics) = c(temp, "lambda_bs_conv_eff")
} else {
    print("Skipping BS-conversion efficiency estimate via lambda alignment ('--with_lambda' not specified)...")
}

######################################################
#  FastQC metrics (pre-trimming)
######################################################

#  Shortened metric names to collect
key_names = c("FQCbasicStats","perBaseQual","perTileQual","perSeqQual",
              "perBaseContent","GCcontent","Ncontent","SeqLengthDist",
              "SeqDuplication","OverrepSeqs","AdapterContent","KmerContent")
                                                                              
#  Get filenames for the FastQC logs (2 per ID for paired-end samples)
if (ncol(manifest) > 3) {
    f = c()
    for (i in 1:length(ids)) {
        f = c(f, 
              paste0(ids[i], '_1_fastqc.log'), 
              paste0(ids[i], '_2_fastqc.log'))
    }
    
    key_names = as.vector(outer(key_names, c('R1', 'R2'), FUN=paste0))
} else {
    f = paste0(ids, '_fastqc.log')
}

if (file.exists(f[1])) {
    print("Adding FastQC metrics (from before any trimming)...")
    
    #  Get the first "column" in each log, and use this ugly code to form it into
    #  a data frame where rows are samples and columns are FastQC metrics
    #  (paired-end samples have one column per read per metric)
    temp_df = data.frame(matrix(unlist(lapply(f, function(x) read.table(x, sep='\t')[,1])),
                                nrow=length(ids),
                                byrow=TRUE))
    colnames(temp_df) = key_names
    
    metrics = cbind(metrics, temp_df)
} else {
    print("Skipping FastQC metrics (no files specified)...")
}

rownames(metrics) = ids
save(metrics, file='metrics.rda')
print("Metrics saved as 'metrics.rda'.")
