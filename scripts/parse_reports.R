library('jaffelab')

#  Get sample IDs
manifest = read.table('arioc_samples.manifest', header = FALSE, stringsAsFactors = FALSE)
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
              
f = paste0(ids, '_alignment.log')
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

filenames = paste0('BME_', ids, '.log')
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
    
filepaths = list.files(pattern='.*\\.f.*q.*_trimming_report\\.txt')

if (length(filepaths) > 0) {
    print('Extracting trimming metrics...')
    
    df_entries = c()
    
    is_first_iter = TRUE
    for (filename in filepaths) {
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
    }
    
    #  Use the correct colnames for paired-end reads
    if (nrow(manifest) > 3) {
        col_names = as.vector(outer(col_names, c('_R1', '_R2'), paste0))
    }
    
    #  Form a data frame
    temp_df = matrix(df_entries, ncol=length(col_names), byrow=TRUE)
    temp_df = as.data.frame(temp_df)
    colnames(temp_df) = col_names
    
    #  Make the integer-containing columns numeric
    numeric_cols = which(!grepl('Sequence', col_names))
    for (i in numeric_cols) {
        temp_df[,i] = as.numeric(temp_df[,i])
    }
    
    metrics = cbind(metrics, temp_df)
} else {
    print("Skipping trimming metrics (no samples appeared to be trimmed)...")
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

rownames(metrics) = ids
save(metrics, file='metrics.rda')
print("Metrics saved as 'metrics.rda'.")
