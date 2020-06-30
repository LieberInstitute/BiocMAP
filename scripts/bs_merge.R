suppressPackageStartupMessages(library('bsseq'))
suppressPackageStartupMessages(library('HDF5Array'))
suppressPackageStartupMessages(library('getopt'))
suppressPackageStartupMessages(library('BiocParallel'))

spec <- matrix(c('dir', 'd', 1, 'character', 'base directory containing objects',
                 'context', 'c', 1, 'character', 'cytosine context of object'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

chrs = readLines(list.files(pattern='chr_names_.*'))
bs_dirs = paste0(opt$dir, '/', chrs)

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
bs_list = list()
for (i in 1:length(chrs)) {
    print(paste0("Loading ", chrs[i], "..."))
    bs_list[[i]] = loadHDF5SummarizedExperiment(bs_dirs[i], prefix=opt$context)
}
gc()

#  Merge entire objects in a delayed 'rbind' call
print("Merging objects (assay merging delayed)...")
bs_big = do.call(rbind, bs_list)
gc()

#  Save the combined result
print("Saving...")
saveHDF5SummarizedExperiment(bs_big, dir=paste0(opt$dir, '/combined'),
                             prefix=opt$context, verbose=TRUE)
print("Done.")

proc.time() - r_time
rm(list=ls())
