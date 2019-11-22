library('jaffelab')

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
              
f = system('ls *_alignment.log', intern=TRUE)
ids = ss(f, '_alignment.log')

row_data = lapply(f, function(filename) {
    log_text = gsub(' ', '', system(paste0('cat ', filename, ' | cut -d " " -f 3-'), intern=TRUE))
    first_row = match("SAMoutput:", log_text) + 1
    log_text = log_text[first_row: (first_row + length(col_names) - 1)]
    
    as.numeric(cleanseVec(ss(log_text,":", 2)))
})

arioc_results = as.data.frame(matrix(unlist(row_data), nrow=length(row_data)))
colnames(arioc_results) = col_names
rownames(ids)

save(arioc_results, file="alignment_results.rda")