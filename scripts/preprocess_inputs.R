#  1. Unzip any gzipped files
#  2. Rename any .fq files to .fastq
#  3. Rename any paired-end suffices to "_1" and "_2"
#  4. Merge any files with same ID

## Required libraries
library('devtools')
library('jaffelab')

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

manifest = read.table('samples.manifest', sep = ' ', header = FALSE, stringsAsFactors = FALSE)

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
for (i in which(exts == ".gz")) {
    command = paste0('gunzip ', filenames[i])
    run_command(command)
    filenames[i] = ss(filenames[i], ".gz")
}

#  Recompute last file extensions
split_name = strsplit(filenames, ".", fixed=TRUE)
exts = sapply(split_name, function(x) x[length(x)])

#  Rename any .fq files to have extension ".fastq" (to simplify file handling
#  in the pipeline)
print("Renaming any '.fq' files to have extension '.fastq'...")
for (i in which(exts == ".fq")) {
    new_name = paste0(ss(filenames[i], '.fq'), '.fastq')
    command = paste('mv', filenames[i], new_name)
    run_command(command)
}

if (paired) {
    print("Replacing any paired suffices with '_1' and '_2'...")
    for (i in 1:nrow(manifest)) {
        new_file1 = paste0(manifest[i, 5], '_1.fastq')
        command = paste('mv', manifest[i, 1], new_file1)
        run_command(command)
        
        new_file2 = paste0(manifest[i, 5], '_2.fastq')
        command = paste('mv', manifest[i, 3], new_file2)
        run_command(command)
    }
}

#  This forms a list, where each element is a vector containing row numbers
#  of the manifest to combine (and each element contains a unique set of rows)
indicesToCombine = list()
for (i in 1:nrow(manifest)) {
    idMatchRows = which(manifest[, ncol(manifest)] == manifest[i, ncol(manifest)])
    if (all(idMatchRows >= i) && length(idMatchRows) > 1) {
        indicesToCombine[[i]] = idMatchRows
    }
}

#  Merge any files that need to be merged (based on having the same ID in the
#  last column of the manifest
print("Merging any files that need to be merged...")
if (paired) {
    for (indices in indicesToCombine) {
        files_to_combine = do.call(paste, as.list(manifest[indices, 1]))
        new_file = paste0(manifest[indices[1], 5], '_R1.fastq')
        command = paste0('cat ', files_to_combine, ' > ', new_file)
        run_command(command)
        command = paste0('rm ', files_to_combine)
        run_command(command)
        
        files_to_combine = do.call(paste, as.list(manifest[indices, 3]))
        new_file = paste0(manifest[indices[1], 5], '_R2.fastq')
        command = paste0('cat ', files_to_combine, ' > ', new_file)
        run_command(command)
        command = paste0('rm ', files_to_combine)
        run_command(command)
    }
} else {
    for (indices in indicesToCombine) {
        files_to_combine = do.call(paste, as.list(manifest[indices, 1]))
        new_file = paste0(manifest[indices[1], 5], '.fastq')
        command = paste0('cat ', files_to_combine, ' > ', new_file)
        run_command(command)
        command = paste0('rm ', files_to_combine)
        run_command(command)
    }
}

#  A sloppy workaround for nextflow not recognizing input files as valid
#  output files (appears to be based on filename- touching the inputs doesn't
#  work as a solution)
print("Giving files temporary new name to ensure placement in channel...")
result_fastqs = system('ls *.fastq', intern=TRUE)
for (f in result_fastqs) {
    file_base = ss(f, ".", fixed=TRUE)
    system(paste0('mv ', f, ' ', file_base, '.fastq_temp'))
}

print("Done with all input preprocessing.")