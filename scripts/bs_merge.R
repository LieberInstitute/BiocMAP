library('bsseq')
library('HDF5Array')
library('getopt')

spec <- matrix(c('dir', 'd', 1, 'character', 'directory to write results to'),
               byrow=TRUE, ncol=5)
opt <- getopt(spec)

bs_objs_CpG = system('ls *_CpG.rda', intern=TRUE)
bs_objs_CpH = system('ls *_CpH.rda', intern=TRUE)
writePath_CpG = paste0(opt$dir, '/bs_combined_CpG.rda')
writePath_CpH = paste0(opt$dir, '/bs_combined_CpH.rda')
writePathAssays = paste0(opt$dir, '/assays.h5')

#  jhpce-specific requirements
BiocParallel::register(MulticoreParam(1))
setAutoBPPARAM(MulticoreParam(1))  # this is necessary, despite above line
h5disableFileLocking()
setAutoBlockSize(1e9)  # not required; slower default is 1e8

###########################################################
#  Merge Bsseq objects containing CpG-context cytosines
###########################################################

print('Combining CpG Bsseq objects...')
load(bs_objs_CpG[1])
BS_full = strandCollapse(BS_CpG[order(ranges(BS_CpG)),])
rm(BS_CpG)
gc()
print("Loaded and strand collapsed first CpG object.")
for (i in 2:length(bs_objs_CpG)) {
  load(bs_objs_CpG[i])
  BS_CpG = strandCollapse(BS_CpG[order(ranges(BS_CpG)),])
  BS_full = rbind(BS_full, BS_CpG)
  print(paste0("Loaded, collapsed, and appended CpG object ", i, "."))
  rm(BS_CpG)
  gc()
}
BSobj = BS_full
rm(BS_full)

print("Writing CpG assays to .h5 file...")
M_h5 = writeHDF5Array(assays(BSobj)$M, writePathAssays, 'M_CpG')
Cov_h5 = writeHDF5Array(assays(BSobj)$Cov, writePathAssays, 'Cov_CpG')

print("Pointing assays to combined CpG object...")
assays(BSobj)$M = M_h5
assays(BSobj)$Cov = Cov_h5

print("Saving remaining CpG object...")
save(BSobj, file=writePath_CpG)

###########################################################
#  Merge Bsseq objects containing CpH-context cytosines
###########################################################

print('Combining CpH Bsseq objects...')
load(bs_objs_CpH[1])
BS_full = BS_CpH
rm(BS_CpH)
gc()
print("Loaded first CpH object.")
for (i in 2:length(bs_objs_CpG)) {
  load(bs_objs_CpH[i])
  BS_full = rbind(BS_full, BS_CpH)
  print(paste0("Loaded and appended CpH object ", i, "."))
  rm(BS_CpH)
  gc()
}
BSobj = BS_full
rm(BS_full)

print("Writing CpH assays to .h5 file...")
M_h5 = writeHDF5Array(assays(BSobj)$M, writePathAssays, 'M_CpH')
Cov_h5 = writeHDF5Array(assays(BSobj)$Cov, writePathAssays, 'Cov_CpH')

print("Pointing assays to combined CpH object...")
assays(BSobj)$M = M_h5
assays(BSobj)$Cov = Cov_h5

print("Saving remaining CpH object...")
save(BSobj, file=writePath_CpH)

rm(list=ls())
print("Done with all tasks.")