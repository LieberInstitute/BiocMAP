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
    
#  String patterns to remove via gsub for every value to put in a data frame:
#  note that patterns are removed in order of occurrence in this list
blacklist = c(" ", "bp", ",", ";Type", "(\\(.*\\))")
cleanseVec <- function(vec) {
    #  Remove blacklisted characters
    for (i in 1:length(blacklist)) {
        vec <- sapply(vec, function(x) gsub(blacklist[i],"",x, perl=TRUE))
    }
    return(vec)
}
    
filepaths = list.files(pattern='.*\\.f.*q.*_trimming_report\\.txt')

if (length(filepaths) > 0) {
    row_starts = c("Total reads processed:",
                   "Reads with adapters:",
                   "Reads written (passing filters):",
                   "Total basepairs processed:",
                   "Quality-trimmed:",
                   "Total written (filtered):",
                   "Sequence:")
    
    clean = function(text_line) {
        text_line = ss(text_line, ":", 2)
        
        #  Remove anything following ";" or "("
        for (flag in c("(", ";")) {
            if (grepl(flag, text_line, fixed=TRUE)) {
                text_line = ss(text_line, flag, fixed=TRUE)
            }
        }
        
        #  Remove these phrases
        for (phrase in c(" ", "bp", ",")) {
            text_line = gsub(phrase, "", text_line)
        }
        
        return(text_line)
    }
    
    parse_report = function(path) {
        report_lines = readLines(path)
        
        return(sapply(row_starts, function(flag) clean(report_lines[grep(flag, report_lines, fixed=TRUE)])))
    }
    
    
    paired = ncol(manifest) > 3
    col_names = gsub(" ", "_", row_starts)
    col_names = sub(":", "", col_names)
    
    if (paired) {
        col_names = as.vector(outer(col_names,
                                    c("_R1", "_R2"),
                                    FUN=paste0))
        trim_metrics = matrix(sapply(filepaths, parse_report), 
                              byrow = TRUE,
                              nrow = length(ids),
                              dimnames = list(ids, col_names))
        
        #  Fix the data type (most columns are integers)
        trim_metrics = as.data.frame(trim_metrics)
        for (i in 1:(ncol(trim_metrics) - 1)) {
            if (i != (ncol(trim_metrics) / 2)) {
            #if (i %% ncol(trim_metrics) < ncol(trim_metrics) - 1) {
                trim_metrics[,i] = as.numeric(trim_metrics[,i])
            }
        }
        
    } else {
        #  Form the matrix of metrics
        trim_metrics = matrix(sapply(filepaths, parse_report), 
                              byrow = TRUE,
                              nrow = length(ids),
                              dimnames = list(ids, col_names))
                              
        #  Fix the data type (most columns are integers)
        trim_metrics = as.data.frame(trim_metrics)
        for (i in 1:(ncol(trim_metrics) - 1)) {
            trim_metrics[,i] = as.numeric(trim_metrics[,i])
        }
    }
    
    metrics = cbind(metrics, trim_metrics)
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
