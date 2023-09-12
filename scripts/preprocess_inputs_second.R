library('jaffelab')

###############################################################################
#  Functions
###############################################################################

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

get_value = function(rules, key, required) {
    this_line = rules[grep(key, ss(rules, '='))]
    if (length(this_line) == 0) {
        if (required) {
            stop(paste0('"', key, '" is a required key in "rules.txt", but was not found.'))
        } else {
            return(NA)
        }
    } else if (length(this_line) > 1) {
        stop(paste0("Key '", key, "' had ", length(this_line), " values in rules.txt."))
    }
    
    return(gsub(' ', '', ss(this_line, '=', 2)))
}

get_paths = function(rules, key, ids, required) {
    value = get_value(rules, key, required)
    
    if (!is.na(value)) {
        pieces = strsplit(value, '[id]', fixed=TRUE)[[1]]
        
        #  Get the paths, substituting in actual ids for each '[id]'
        paths = ''
        for (i in 1:(length(pieces)-1)) {
            paths = paste0(paths, pieces[i], ids)
        }
        paths = paste0(paths, pieces[length(pieces)])
        
        return(paths)
    } else {
        return(NA)
    }
}

rename_links = function(glob, orig_paths_long, sym_paths_short, end_name, key, expected_matches, missing_ok=FALSE) {
    total_ids_matched = 0
    
    for (i in 1:length(glob)) {
        matches = regexpr(pattern = glob[i], orig_paths_long)
        num_matches = length(which(matches == 1))
        
        if (missing_ok && (num_matches == 0)) {
            next
        } else if (num_matches != expected_matches) {
            stop(paste0(num_matches, ' matches found for "', key , '" per ID when ', expected_matches, ' matches were expected.'))
        }
        
        total_ids_matched = total_ids_matched + 1
        
        sym_path = sym_paths_short[matches == 1]
        
        if (length(sym_path) == 1) {
            run_command(paste('mv', sym_path, paste0(ids[i], end_name)))
        } else {
            #  Assumes 'list.files' will return mates in order!
            run_command(paste('mv', sym_path[1], paste0(ids[i], '_1', end_name)))
            run_command(paste('mv', sym_path[2], paste0(ids[i], '_2', end_name)))
        }
    }
    
    if (total_ids_matched == 0) {
        stop("Didn't find any files for '", key , "' (expected at least 1 ID to have a matching file).")
    }
}

get_ext = function(filename) {
    num_fields = length(strsplit(filename, '\\.')[[1]])
    
    last_field = ss(filename, '\\.', num_fields)
    if (last_field == 'bai') {
        ext = paste0(ss(filename, '\\.', num_fields - 1), '.', last_field)
        if (ext != 'bam.bai') {
            stop(paste0('Found an unknown file type: ".', ext, '".'))
        }
        return (ext)
    } else {
        return (last_field)
    }
}

###############################################################################
#  Rename logs with meaningful, unique filenames
###############################################################################

#  Read in 'rules.txt'
rules = readLines('rules.txt')

if (length(grep('#', rules)) > 0) {
    rules = rules[-grep('#', rules)]
}

#  Read in 'samples.manifest'
manifest = read.table(
    'samples.manifest', header = FALSE, stringsAsFactors = FALSE
)
ids = manifest[,ncol(manifest)]
paired = ncol(manifest) > 3

#  Name of the symbolic links present in the working directory, full path to
#  the files each link points to, and basename of these files, respectively
sym_paths_short = list.files(pattern = '.*\\.$')
orig_paths_long = Sys.readlink(sym_paths_short)
orig_paths_short = basename(orig_paths_long)

#  For each type of log:
#    1. Grab the glob specified in rules.txt for the key
#    2. See which files in the WD match the glob
#    3. Rename those files appropriately

#  Arioc logs
print('Looking for Arioc logs...')
glob = get_paths(rules, 'arioc_log', ids, TRUE) # now we have one glob per ID
rename_links(glob, orig_paths_long, sym_paths_short, '_arioc.log', 'arioc_log', 1)

#  XMC logs
glob = get_paths(rules, 'xmc_log', ids, FALSE) # now we have one glob per ID
if (!is.na(glob[1])) {
    print('Looking for XMC logs, since the "xmc_log" key was specified...')
    rename_links(glob, orig_paths_long, sym_paths_short, '_xmc.log', 'xmc_log', 1)
}

#  FastQC logs
print('Looking for FastQC logs...')
glob_last = get_paths(rules, 'fastqc_log_last', ids, FALSE)
glob_first = get_paths(rules, 'fastqc_log_first', ids, FALSE)

num_matches = ifelse(paired, 2, 1)
if (is.na(glob_first[1]) && !is.na(glob_last[1])) {
    #  If only 'last' FastQC logs exist, these will be the ones we parse
    rename_links(glob_last, orig_paths_long, sym_paths_short, '_fastqc.log', 'fastqc_log_last', num_matches)
} else if (!is.na(glob_last[1])) {
    #  For simplicity given the existing code/functions, all "first fastqc logs"
    #  are renamed, and renaming for any "last fastqc logs" is performed afterward,
    #  thus overwriting "first fastqc logs" for applicable IDs
    rename_links(glob_first, orig_paths_long, sym_paths_short, '_fastqc.log', 'fastqc_log_first', num_matches)
    rename_links(glob_last, orig_paths_long, sym_paths_short, '_fastqc.log', 'fastqc_log_last', num_matches, TRUE)
}

#  Trim Galore reports
glob = get_paths(rules, 'trim_report', ids, FALSE) # now we have one glob per ID
if (!is.na(glob[1])) {
    print('Looking for Trim Galore logs, since the "trim_report" key was specified...')
    rename_links(glob, orig_paths_long, sym_paths_short, '_trim_report.log', 'trim_report', 1, TRUE)
}

#  SAM or (BAM and .bam.bai) files
print('Looking for alignment-related files...')
glob = get_paths(rules, 'sam', ids, FALSE) # now we have one glob per ID

temp = sapply(orig_paths_short, get_ext) == 'sam'
orig_sams = orig_paths_long[temp]

#  Note that 'rename_links' would catch the problematic case where
#  0 < length(orig_sams) < length(ids)
if (length(orig_sams) > 0) {
    print('Found at least one SAM file; looking for the remaining ones...')
    sym_sams = sym_paths_short[temp]
    rename_links(glob, orig_sams, sym_sams, '.sam', 'sam', 1)
} else {
    print("Didn't find any SAM files; looking for BAMs instead...")
    temp = sapply(orig_paths_short, get_ext) == 'bam'
    rename_links(glob, orig_paths_long[temp], sym_paths_short[temp], '.bam', 'sam', 1)
    
    print('Looking for indices for these BAMs (".bam.bai" files)...')
    temp = sapply(orig_paths_short, get_ext) == 'bam.bai'
    rename_links(glob, orig_paths_long[temp], sym_paths_short[temp], '.bam.bai', 'sam', 1)
}

###############################################################################
#  Verify the files in the manifest have known extensions
###############################################################################

#  Nextflow passed files from the manifest to the working directory by this
#  point, so relative paths should be used (and in fact must be used when
#  running with docker- absolute paths are not mounted in this container)
manifest[, 1] <- basename(manifest[, 1])
if (paired) manifest[, 3] <- basename(manifest[, 3])

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
    stop("Unrecognized FASTQ filename extension (should be 'fastq.gz', 'fq.gz', 'fastq' or 'fq').")
}

if (length(unique(actual_exts)) > 1) {
    stop("We require consistent file extensions for FASTQ files among samples and between reads. The extensions 'fastq.gz', 'fq.gz', 'fastq' or 'fq' are all accepted, but all files must use the same extension.")
}

###############################################################################
#  Perform merging and renaming of files for use in the pipeline
###############################################################################

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

#   Verify file extensions match for different files (rows) with the same sample
#   ID
error_message <- "FASTQ file extensions must match among all files for a single sample."
if (!all(sapply(indicesToCombine, function(x) length(unique(actual_exts[x])) == 1))) {
    stop(error_message)
}
if (paired && !all(sapply(indicesToCombine, function(x) length(unique(actual_exts[nrow(manifest) + x])) == 1))) {
    stop(error_message)
}

if (paired) {
    #  Merge files that require merging
    print("Merging any files that need to be merged...")
    for (indices in indicesToCombine) {
        #  Determine file extension for the merged file to have
        first_ext <- actual_exts[indices[1]]

        #  Do the file merging
        files_to_combine <- do.call(paste, as.list(manifest[indices, 1]))
        new_file <- paste0(manifest[indices[1], 5], "_1.", first_ext)
        command <- paste("cat", files_to_combine, ">", new_file)
        run_command(command)

        #  Note that 'basename' was called on the manifest, so unintended
        #  application of this command will not do harm outside the nextflow
        #  temporary directory
        command <- paste("rm", files_to_combine)
        run_command(command)

        files_to_combine <- do.call(paste, as.list(manifest[indices, 3]))
        new_file <- paste0(manifest[indices[1], 5], "_2.", first_ext)
        command <- paste("cat", files_to_combine, ">", new_file)
        run_command(command)

        #  Note that 'basename' was called on the manifest, so unintended
        #  application of this command will not do harm outside the nextflow
        #  temporary directory
        command <- paste("rm", files_to_combine)
        run_command(command)
    }

    #  Rename any remaining files by their associated sampleID, using the
    #  paired suffices _1 and _2.
    print("Renaming remaining files for handling in the pipeline...")
    if (length(unlist(indicesToCombine)) > 0) {
        remaining_rows <- (1:nrow(manifest))[-unlist(indicesToCombine)]
    } else {
        remaining_rows <- 1:nrow(manifest)
    }
    for (index in remaining_rows) {
        first_ext <- actual_exts[index]

        new_file <- paste0(manifest[index, 5], "_1.", first_ext)
        if (manifest[index, 1] != new_file) {
            command <- paste("mv", manifest[index, 1], new_file)
            run_command(command)
        }

        new_file <- paste0(manifest[index, 5], "_2.", first_ext)
        if (manifest[index, 3] != new_file) {
            command <- paste("mv", manifest[index, 3], new_file)
            run_command(command)
        }
    }
} else {
    print("Merging any files that need to be merged...")
    for (indices in indicesToCombine) {
        #  Do the file merging
        files_to_combine <- do.call(paste, as.list(manifest[indices, 1]))
        new_file <- paste0(manifest[indices[1], 3], ".", actual_exts[indices[1]])
        command <- paste("cat", files_to_combine, ">", new_file)
        run_command(command)

        #  Note that 'basename' was called on the manifest, so unintended
        #  application of this command will not do harm outside the nextflow
        #  temporary directory
        command <- paste("rm", files_to_combine)
        run_command(command)
    }

    #  Rename remaining files by their associated sampleID
    print("Renaming remaining files for handling in the pipeline...")
    if (length(unlist(indicesToCombine)) > 0) {
        remaining_rows <- (1:nrow(manifest))[-unlist(indicesToCombine)]
    } else {
        remaining_rows <- 1:nrow(manifest)
    }
    for (index in remaining_rows) {
        new_file <- paste0(manifest[index, 3], ".", actual_exts[index])
        if (manifest[index, 1] != new_file) {
            command <- paste("mv", manifest[index, 1], new_file)
            run_command(command)
        }
    }
}


print("Done all tasks.")
