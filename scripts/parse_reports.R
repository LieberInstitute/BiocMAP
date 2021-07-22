library('jaffelab')

get_value = function(rules, key, required) {
    this_line = rules[grep(key, ss(rules, '='))]
    
    return(gsub(' ', '', ss(this_line, '=', 2)))
}

#  Get sample IDs
rules = readLines('rules.txt')

if (length(grep('#', rules)) > 0) {
    rules = rules[-grep('#', rules)]
}

manifest = read.table(get_value(rules, 'manifest'), header = FALSE,
                      stringsAsFactors = FALSE)
ids = manifest[,ncol(manifest)]
paired = ncol(manifest) > 3

######################################################
#  Arioc
######################################################

print("Extracting Arioc metrics...")

f = paste0(ids, '_arioc.log')
stopifnot(all(file.exists(f)))

parse_rows = function(log_text, blacklist) {
    #  Get the "value" portion of each row
    log_text = ss(log_text, ':', 2)
    
    #  Remove blacklisted characters
    for (b in blacklist) {
        log_text = gsub(b, "", log_text)
    }
    
    #  Return values in each row
    return(as.numeric(log_text))
}

#  Return a numeric vector of values given a single path to an alignment log
parse_log = function(filename, expected_keys) {
    #  Process log into a "machine-friendly" format with the relevant information
    log_text = system(paste0('cat ', filename, ' | cut -d " " -f 3-'), intern=TRUE)
    log_text = gsub(' ', '', log_text)
    first_row = match("SAMoutput:", log_text) + 1
    log_text = log_text[first_row: (first_row + length(expected_keys) - 1)]
    
    #  Verify the information is as expected
    actual_keys = ss(log_text, ":")
    stopifnot(all(expected_keys == actual_keys))
    
    if (paired) {
        #  Manually parse TLEN-related metrics
        manual_row = match("TLENmode(mean,sd,skewness)", actual_keys)
        tlen_stats = ss(log_text[manual_row], ":", 2)
        tlen_mode = ss(tlen_stats, '\\(', 1)
        
        tlen_stats = gsub('\\)', '', ss(tlen_stats, '\\(', 2))
        tlen_mean = as.numeric(ss(tlen_stats, ',', 1))
        tlen_sd = as.numeric(ss(tlen_stats, ',', 2))
        tlen_skewness = as.numeric(ss(tlen_stats, ',', 3))
        log_text = log_text[-manual_row]
        actual_keys = actual_keys[-manual_row]
    }
    
    #  Manually parse error rate
    manual_row = match("errorrate(primarymappings)", actual_keys)
    err_stats = ss(log_text[manual_row], ":", 2)
    err_rate = as.numeric(gsub("%", "", err_stats)) / 100
    log_text = log_text[-manual_row]
    actual_keys = actual_keys[-manual_row]
    
    #  Parse all other keys
    blacklist = c("\\(.*\\)")
    values = parse_rows(log_text, blacklist)
    
    if (paired) {
        values = c(values, tlen_mean, tlen_sd, tlen_skewness, err_rate)
    } else {
        values = c(values, err_rate)
    }
    
    return(values)
}

#  Define the expected literal keys to observe in the log ("expected_keys")
#  and the name of the column we wish to associate with those keys
#  (note some "special" columns are added in the "parse_log" function
if (paired) {
    col_names = c(
        "pairs", "conc_pairs_total", "conc_pairs_1_mapping",
        "conc_pairs_many_mappings", "disc_pairs", "unmapped_total_pairs",
        "unmapped_rejected_pairs", "unmapped_other_pairs",
        "mates_not_in_paired_maps_total", "mates_NIPM_with_no_maps",
        "mates_NIPM_with_1_map", "mates_NIPM_with_many_maps",
        "total_mapped_mates", "duplicate_maps", "maxQlen",
        "max_diag_band_width", "TLEN_disc_pairs", "TLEN_mean", "TLEN_sd",
        "TLEN_skewness", "err_rate_primary_maps"
    )
    
    expected_keys = c(
        "pairs", "concordantpairs", "with1mapping", "with2ormoremappings",
        "discordantpairs", "unmappedpairs", "twomatesmapped(\"rejected\")",
        "0or1matesmapped", "matesnotinpairedmappings", "withnomappings",
        "with1mapping", "with2ormoremappings", "totalmappedmates",
        "duplicatemappings(unreported)", "maximumQlength",
        "maximumdiagonalbandwidth", "TLENmode(mean,sd,skewness)",
        "TLEN-discordantpairs", "errorrate(primarymappings)"
    )
} else {
    col_names = c(
        "reads", "mapped_reads_total", "mapped_reads_1_map",
        "mapped_reads_many_maps", "unmapped_reads", "duplicate_maps",
        "err_rate_primary_maps"
    )
    
    expected_keys = c(
        "reads", "mappedreads", "with1mapping", "with2ormoremappings",
        "unmappedreads", "duplicatemappings(unreported)",
        "errorrate(primarymappings)"
    )
}

metrics = as.data.frame(
    matrix(
        unlist(lapply(f, parse_log, expected_keys)),
        ncol=length(col_names),
        byrow=TRUE,
        dimnames=list(ids, col_names)
    )
)

######################################################
#  MethylDackel or BME logs
######################################################

print('Extracting methylation extraction metrics...')
filenames = paste0('methyl_extraction_', ids, '.log')
stopifnot(all(file.exists(filenames)))

percs = lapply(filenames, function(f) {
    raw_chars = system(paste0('cat ', f, ' | grep -E "C methylated in C[pH]?. context:" | cut -d ":" -f 2'), intern=TRUE)
    return(as.numeric(sub("[:space:]|%", "", raw_chars)))
})
mat = matrix(unlist(percs), dimnames = list(ids, c("perc_M_CpG", "perc_M_CHG", "perc_M_CHH")), byrow=TRUE, ncol=3)

metrics = cbind(metrics, mat)

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
#  FastQC metrics (after trimming where applicable)
######################################################

#  Possible columns we expect in FastQC summaries
fastqc_stat_names = c(
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
    "Kmer Content"
)

#  Parse a vector of FastQC summaries, containing each sample in the order
#  present in 'ids'. Return a matrix of statistics
parse_summary <- function(summary_paths, ids) {
    temp_rows <- lapply(summary_paths, function(x) {
        temp <- readLines(x)
        vals <- ss(temp, "\t")
        names(vals) <- ss(temp, "\t", 2)

        stopifnot(all(names(vals) %in% fastqc_stat_names))
        vals <- vals[fastqc_stat_names]
    })

    actual_colnames <- gsub(" ", "_", tolower(fastqc_stat_names))
    stat_mat <- matrix(unlist(temp_rows),
        ncol = length(fastqc_stat_names),
        byrow = TRUE,
        dimnames = list(ids, actual_colnames)
    )

    return(stat_mat)
}

if (file.exists(paste0(ids[1], '_1_fastqc.log')) || file.exists(paste0(ids[1], '_fastqc.log'))){
    print("Adding FastQC metrics (from before any trimming)...")
    
    if (paired) {
        summary_paths_1 = paste0(ids, '_1_fastqc.log')
        summary_paths_2 = paste0(ids, '_2_fastqc.log')
        
        #  Get stat values for each mate
        stat_mat <- parse_summary(summary_paths_1, ids)
        stat_mat_2 <- parse_summary(summary_paths_2, ids)
    
        #  Combine into a single matrix and simplify some values
        stat_mat[] <- paste(stat_mat, stat_mat_2, sep = "/")
        stat_mat[stat_mat == "PASS/PASS"] <- "PASS"
        stat_mat[stat_mat == "WARN/WARN"] <- "WARN"
        stat_mat[stat_mat == "FAIL/FAIL"] <- "FAIL"
        stat_mat[stat_mat == "NA/NA"] <- "NA"
    
        metrics <- cbind(metrics, stat_mat)
    } else {
        f = paste0(ids, '_fastqc.log')
        metrics <- cbind(metrics, parse_summary(f, ids))
    }
} else {
    print("Skipping FastQC metrics (no files specified)...")
}

rownames(metrics) = ids
save(metrics, file='metrics.rda')
print("Metrics saved as 'metrics.rda'.")
