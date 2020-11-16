library('jaffelab')

###################################################
#  Functions
###################################################

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

get_value = function(rules, key) {
    this_line = rules[grep(key, ss(rules, '='))]
    if (length(this_line) != 1) {
        stop(paste0("Key '", key, "' had ", length(this_line), " values in rules.txt."))
    }
    
    return(gsub(' ', '', ss(this_line, '=', 2)))
}

form_links = function(rules, key, ids, end_name) {
    value = get_value(rules, key)
    if (value != "NA") {
        pieces = strsplit(value, '[id]', fixed=TRUE)[[1]]
        
        #  Get the paths, substituting in actual ids for each '[id]'
        paths = ''
        for (i in 1:(length(pieces)-1)) {
            paths = paste0(paths, pieces[i], ids)
        }
        paths = paste0(paths, pieces[length(pieces)])
        stopifnot(length(paths) == length(ids))
        
        #  Symbolically link each file into the current working directory
        for (i in 1:length(ids)) {
            command = paste0('ln -s ', paths[i], ' ', ids[i], end_name)
            run_command(command)
        }
    }
}

###################################################
#  Link logs into the working directory
###################################################

rules = readLines('rules.txt')
rules = rules[-grep('#', rules)]

manifest = read.table(get_value(rules, 'manifest'), header = FALSE,
                      stringsAsFactors = FALSE)
ids = manifest[,ncol(manifest)]
paired = ncol(manifest) > 3


#  Arioc SAMs and logs
form_links(rules, 'sam', ids, '.sam')
form_links(rules, 'arioc_log', ids, '_arioc.log')

#  XMC logs
form_links(rules, 'xmc_log', ids, '_xmc.log')

#  BME logs
form_links(rules, 'bme_log', ids, '_bme.log')

#  Trim Galore reports
form_links(rules, 'trim_report_r1', ids, '_trim_report_r1.txt')
form_links(rules, 'trim_report_r2', ids, '_trim_report_r2.txt')


############################################################
#  Verify the FASTQs in the manifest have known extensions
############################################################

print("Verifying file extensions from the manifest are valid...")

#  Get the FASTQ filenames, as declared in samples.manifest
filenames = manifest[,1]
if (paired) filenames = c(filenames, manifest[,3])

is_zipped = grepl(".gz", filenames, fixed=TRUE)

valid_exts = c('fastq.gz', 'fq.gz', 'fastq', 'fq')
actual_exts = sapply(1:length(filenames), function(i) {
    temp = strsplit(sub(".gz", "", filenames[i]), ".", fixed=TRUE)[[1]]
    if (is_zipped[i]) {
        return(paste0(temp[length(temp)], '.gz'))
    } else {
        return(temp[length(temp)])
    }
})

if (!all(actual_exts %in% valid_exts)) {
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}

if (paired && any(actual_exts[1:nrow(manifest)] != actual_exts[(nrow(manifest)+1):length(actual_exts)])) {
    stop("A given pair of reads must have the same file extensions.")
}


####################################################################
#  Perform merging and renaming of FASTQs for use in the pipeline
####################################################################

#  This forms a list, where each element is a vector containing row numbers
#  of the manifest to combine (and each element contains a unique set of rows)
indicesToCombine = list()
for (i in 1:nrow(manifest)) {
    idMatchRows = which(manifest[, ncol(manifest)] == manifest[i, ncol(manifest)])
    if (all(idMatchRows >= i) && length(idMatchRows) > 1) {
        indicesToCombine[[i]] = idMatchRows
    }
}
indicesToCombine = indicesToCombine[!sapply(indicesToCombine, is.null)]

if (paired) {
    #  Merge files that require merging
    print("Merging any files that need to be merged...")
    for (indices in indicesToCombine) {
        #  Determine file extension for the merged file to have
        first_ext = actual_exts[indices[1]]
        
        #  Do the file merging
        files_to_combine = do.call(paste, as.list(manifest[indices, 1]))
        new_file = paste0(manifest[indices[1], 5], '_1.', first_ext)
        command = paste('cat', files_to_combine, '>', new_file)
        run_command(command)
        
        files_to_combine = do.call(paste, as.list(manifest[indices, 3]))
        new_file = paste0(manifest[indices[1], 5], '_2.', first_ext)
        command = paste('cat', files_to_combine, '>', new_file)
        run_command(command)
    }
    
    #  Symbolically link any remaining files: this renames the files
    #  by their associated sampleID, and uses the paired suffices _1 and _2.
    print("Renaming and symbolically linking files for handling in the pipeline...")
    if (length(unlist(indicesToCombine)) > 0) {
        remaining_rows = (1:nrow(manifest))[-unlist(indicesToCombine)]
    } else {
        remaining_rows = 1:nrow(manifest)
    }
    for (index in remaining_rows) {
        first_ext = actual_exts[index]
        
        new_file = paste0(manifest[index, 5], '_1.', first_ext)
        command = paste('ln -s', manifest[index, 1], new_file)
        run_command(command)
        
        new_file = paste0(manifest[index, 5], '_2.', first_ext)
        command = paste('ln -s', manifest[index, 3], new_file)
        run_command(command)
    }

} else {
    print("Merging any files that need to be merged...")
    for (indices in indicesToCombine) {
        #  Do the file merging
        files_to_combine = do.call(paste, as.list(manifest[indices, 1]))
        new_file = paste0(manifest[indices[1], 3], '.', actual_exts[indices[1]])
        command = paste('cat', files_to_combine, '>', new_file)
        run_command(command)
    }
    
    #  Symbolically link any remaining files; also renames the files
    #  by their associated sampleID
    print("Renaming and symbolically linking files for handling in the pipeline...")
    if (length(unlist(indicesToCombine)) > 0) {
        remaining_rows = (1:nrow(manifest))[-unlist(indicesToCombine)]
    } else {
        remaining_rows = 1:nrow(manifest)
    }
    for (index in remaining_rows) {
        new_file = paste0(manifest[index, 3], '.', actual_exts[index])
        command = paste('ln -s', manifest[index, 1], new_file)
        run_command(command)
    }
}

print("Done all tasks.")
