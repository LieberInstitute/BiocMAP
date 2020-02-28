## Required libraries
library('devtools')
library('jaffelab')

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

manifest = read.table('samples.manifest', header = FALSE, stringsAsFactors = FALSE)

## Is the data paired end?
paired = ncol(manifest) > 3

#  Get the FASTQ filenames, as declared in samples.manifest
filenames = manifest[,1]
if (paired) filenames = c(filenames, manifest[,3])

#  Get the last file extension (eg. ".gz" for "a.fastq.gz")
split_name = strsplit(filenames, ".", fixed=TRUE)
exts = sapply(split_name, function(x) x[length(x)])

#  Decompress any gzipped files and adjust their filenames to not have
#  the ".gz" ending (Arioc needs unzipped files)
print("Unzipping any gzipped reads...")
for (i in which(exts == "gz")) {
    #  Remove ".gz" and make sure suffix is ".fastq"
    new_name = sub("\\.fq", "\\.fastq", basename(ss(filenames[i], "\\.gz")))
    
    command = paste0('gunzip -c ', filenames[i], ' > ', new_name)
    run_command(command)
    filenames[i] = new_name
}

print("Symbolically linking remaining FASTQs into the working directory...")
for (i in which(exts != "gz")) {
    stopifnot(exts[i] %in% c("fq", "fastq"))
    new_name = sub("\\.fq", "\\.fastq", basename(filenames[i]))
    
    command = paste('ln -s', filenames[i], new_name)
    run_command(command)
    filenames[i] = new_name
}

#  Recompute last file extensions
split_name = strsplit(filenames, ".", fixed=TRUE)
exts = sapply(split_name, function(x) x[length(x)])


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

if (length(unlist(indicesToCombine)) > 0) {
    remaining_rows = (1:nrow(manifest))[-unlist(indicesToCombine)]
} else {
    remaining_rows = 1:nrow(manifest)
}

#  Merge any files that need to be merged (based on having the same ID in the
#  last column of the manifest
print("Merging any files that need to be merged...")
if (paired) {
    for (indices in indicesToCombine) {
        files_to_combine = do.call(paste, as.list(filenames[indices]))
        new_file = paste0(manifest[indices[1], 5], '_1.fastq')
        command = paste0('cat ', files_to_combine, ' > ', new_file)
        run_command(command)
        command = paste0('rm ', files_to_combine)
        run_command(command)
        
        files_to_combine = do.call(paste, as.list(filenames[indices + nrow(manifest)]))
        new_file = paste0(manifest[indices[1], 5], '_2.fastq')
        command = paste0('cat ', files_to_combine, ' > ', new_file)
        run_command(command)
        command = paste0('rm ', files_to_combine)
        run_command(command)
    }
    
    print("Renaming files for handling in the pipeline...")
    for (index in remaining_rows) {
        new_file = paste0(manifest[index, 5], '_1.fastq')
        command = paste('mv', filenames[index], new_file)
        run_command(command)
        
        new_file = paste0(manifest[index, 5], '_2.fastq')
        command = paste('mv', filenames[index + nrow(manifest)], new_file)
        run_command(command)
    }
    
    #  Rewrite a manifest to reflect the file name changes and any merging
    print("Constructing a new manifest to reflect these changes...")
    first_reads = system('ls *_1.f*q*', intern=TRUE)
    ids = ss(first_reads, '_1\\.')
    new_man = paste(first_reads,
                    0,
                    system('ls *_2.f*q*', intern=TRUE),
                    0,
                    ids)
    writeLines(new_man, con="arioc_samples.manifest")
} else {
    for (indices in indicesToCombine) {
        files_to_combine = do.call(paste, as.list(filenames[indices]))
        new_file = paste0(manifest[indices[1], 3], '.fastq')
        command = paste0('cat ', files_to_combine, ' > ', new_file)
        run_command(command)
        command = paste0('rm ', files_to_combine)
        run_command(command)
    }
    
    print("Renaming files for handling in the pipeline...")
    for (index in remaining_rows) {
        new_file = paste0(manifest[index, 3], '.fastq')
        command = paste('mv', filenames[index], new_file)
        run_command(command)
    }
    
    #  Rewrite a manifest to reflect the file name changes and any merging
    print("Constructing a new manifest to reflect these changes...")
    reads = basename(system('ls *.f*q*', intern=TRUE))
    ids = ss(reads, '.', fixed=TRUE)
    new_man = paste(reads, 0, ids)
    writeLines(new_man, con="arioc_samples.manifest")
}
