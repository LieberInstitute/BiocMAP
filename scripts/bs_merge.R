suppressPackageStartupMessages(library('bsseq'))
suppressPackageStartupMessages(library('HDF5Array'))
suppressPackageStartupMessages(library('getopt'))
suppressPackageStartupMessages(library('BiocParallel'))

spec <- matrix(
    c(
        'in_dir', 'i', 1, 'character', 'base directory containing input objects',
        'out_dir', 'o', 1, 'character', 'base directory for output objects',
        'context', 'c', 1, 'character', 'cytosine context of object'
    ),
    byrow = TRUE,
    ncol = 5
)
opt <- getopt(spec)

chrs = readLines(list.files(pattern='chr_names_.*'))
in_dirs = file.path(opt$in_dir, chrs, opt$context)

#  jhpce-specific requirements
BiocParallel::register(MulticoreParam(1))
setAutoBPPARAM(MulticoreParam(1))  # this is necessary, despite above line
h5disableFileLocking()
setAutoBlockSize(1e9)  # not required; slower default is 1e8

###########################################################
#  Main
###########################################################

r_time = proc.time()

#  Load objects (minus their assays) into memory
print(paste("Beginning merging for", opt$context, "context..."))
    
bs_list = list()
for (i in 1:length(chrs)) {
    print(paste0("Loading ", chrs[i], "..."))
    bs_list[[i]] = loadHDF5SummarizedExperiment(in_dirs[i])
}
gc()

#  Merge entire objects in a delayed 'rbind' call
print("Merging objects (assay merging delayed)...")
bs_big = do.call(rbind, bs_list)
gc()

#  Add colData
print("Adding colData...")
load('metrics.rda')
metrics = metrics[colnames(bs_big), ]
colData(bs_big) = DataFrame(metrics)

#  Save the combined result
print("Saving...")
saveHDF5SummarizedExperiment(
    bs_big, dir=opt$out_dir, prefix=opt$context, verbose=TRUE, replace=TRUE
)
print("Done.")

proc.time() - r_time
rm(list=ls())
