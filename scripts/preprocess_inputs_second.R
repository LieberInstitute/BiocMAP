library('jaffelab')

###################################################
#  Functions
###################################################

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

#  Given an absolute file path (character vector of length 1) as a glob
#  expression, return a vector of corresponding filepaths.
process_glob = function(path_glob, key) {
    #  Expand glob and check that the correct number of files were matched
    paths = Sys.glob(path_glob)
    if (length(paths) == 0) {
        stop(paste0('Paths for key "',
                   key,
                   '" in "rules.txt" matched no files.'))
    } else if (key == 'sam') {
        suff = substr(basename(paths), 
                      nchar(basename(paths)) - 3,
                      nchar(basename(paths)))
        using_bams = length(paths) == 2 && all(suff %in% c('.bam', '.bai'))
        using_sams = length(paths) == 1 && all(suff == '.sam')
        if (!using_bams && !using_sams) {
            stop("Incorrect path specification for key 'sam'. Either specify a glob matching '.bam' and '.bam.bai' files, or just match '.sam' files.")
        }
    } else if (key != 'fastqc_log_last' && key != 'fastqc_log_first' && length(paths) > 1) {
        stop(paste0('Paths for key "',
                   key,
                   '" in "rules.txt" matched more than 1 file.'))
    } else if (key != 'fastqc_log_last' && key != 'fastqc_log_first' && !(length(paths) %in% c(1, 2, 4))) {
        stop(paste0('Paths for key "',
                   key,
                   '" in "rules.txt" matched an invalid number of files: ',
                   length(paths)))
    }
    
    return(paths)
}

get_paths = function(rules, key, ids, required) {
    value = get_value(rules, key, required)
    
    if (!is.na(value)) {
        pieces = strsplit(value, '[id]', fixed=TRUE)[[1]]
        
        #  Get the paths, substituting in actual ids for each '[id]', and
        #  evaluating any glob syntax
        
        paths = ''
        for (i in 1:(length(pieces)-1)) {
            paths = paste0(paths, pieces[i], ids)
        }
        paths = paste0(paths, pieces[length(pieces)])
        paths = unlist(lapply(paths, process_glob, key))
        
        return(paths)
    } else {
        return(NA)
    }
}

#  Given absolute paths to the logs for a particular key in 'rules.txt', stop
#  if any files don't exist or if the wrong number of paths exist
check_paths = function(paths, key, ids) {
    if (!is.na(paths)) {
        #  Verify legitimacy of paths
        if (length(paths) != length(ids) && length(paths) != 2 * length(ids)) {
            stop('Improper number of total logs found for key "', key, '".')
        }
        if (! all(file.exists(paths))) {
            stop(paste0('Some or all files for key "', key, '" do not exist.'))
        }
    }
}

form_links = function(paths, ids, end_name) {
    #  Whether 2 logs exist per ID
    paired = length(paths) == 2 * length(ids)
    
    #  Symbolically link each file into the current working directory
    for (i in 1:length(ids)) {
        if (paired) {
            command = paste0('ln -s ', paths[2*i - 1], ' ', ids[i], '_1', end_name)
            run_command(command)
            
            command = paste0('ln -s ', paths[2*i], ' ', ids[i], '_2', end_name)
            run_command(command)
        } else {
            command = paste0('ln -s ', paths[i], ' ', ids[i], end_name)
            run_command(command)
        }
    }
}

process_key = function(rules, key, ids, end_name, required) {
    paths = get_paths(rules, key, ids, required)
    check_paths(paths, key, ids)
    form_links(paths, ids, end_name)
}

###################################################
#  Link logs into the working directory
###################################################

rules = readLines('rules.txt')

if (length(grep('#', rules)) > 0) {
    rules = rules[-grep('#', rules)]
}

manifest = read.table(get_value(rules, 'manifest', TRUE), header = FALSE,
                      stringsAsFactors = FALSE)
ids = manifest[,ncol(manifest)]
paired = ncol(manifest) > 3


#  Arioc logs
process_key(rules, 'arioc_log', ids, '_arioc.log', TRUE)

#  XMC logs
process_key(rules, 'xmc_log', ids, '_xmc.log', FALSE)

#  Trim Galore reports
process_key(rules, 'trim_report', ids, '_trim_report.log', TRUE)

#  FastQC logs
last_paths = get_paths(rules, 'fastqc_log_last', ids, TRUE)
first_paths = get_paths(rules, 'fastqc_log_first', ids, FALSE)

if (is.na(first_paths)) {
    combined_paths = last_paths
    check_paths(combined_paths, 'fastqc_log_last', ids)
} else {
    combined_paths = c(first_paths[!(first_paths %in% last_paths)], last_paths)
    check_paths(paths, 'combined fastqc_log_first and fastqc_log_last', ids)
}

form_links(combined_paths, ids, '_fastqc.log')

#  SAM or (BAM and .bam.bai) files
sam_paths = get_paths(rules, 'sam', ids, TRUE)
check_paths(sam_paths, 'sam', ids)
if (length(sam_paths) == 2 * length(ids)) {
    suff = substr(sam_paths, nchar(sam_paths) - 3, nchar(sam_paths))
    
    #  Note that ids and files to link are guaranteed to be lined up correctly
    form_links(sam_paths[suff == '.bam'], ids, '.bam')
    form_links(sam_paths[suff == '.bai'], ids, '.bam.bai')
} else {
    form_links(sam_paths, ids, '.sam')
}

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