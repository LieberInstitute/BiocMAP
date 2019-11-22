#  Read cytosine reports and produce 2 HDF5-backed BSseq objects (CpG and nonCpG context).

library('bsseq')
library('GenomicRanges')
library('data.table')
library('getopt')
library('HDF5Array')
library('jaffelab')

################################################################
#  Convenience variables
################################################################

spec <- matrix(c(
    'chr', 's', 1, 'character', 'chromosome to subset',
    'region', 'r', 1, 'character', 'brain region'
    'cores', 'c', 1, 'integer', 'number of cores to utilize'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

assayPath = paste0('assays_', opt$chr, '.h5')
tempDir = paste0('temp_hd5/', opt$chr, '/')

#  Ensure clusters like JHPCE don't improperly allocate resources
BiocParallel::register(MulticoreParam(1))
setAutoBPPARAM(MulticoreParam(1))  # this is necessary, despite above line
h5disableFileLocking()
data.table::setDTthreads(opt$cores)
setAutoBlockSize(1e9)  # for speed- default is 1e8

################################################################
#  Construct one BSseq object for each chromosome (later to be
#  split by CpG/nonCpG context). 
################################################################

reportFiles = system('ls *.CX_report.txt', intern=TRUE)
ids = ss(reportFiles, '.', fixed=TRUE)

#  Read the reports
print('Beginning read.bismark...')
BSobj = read.bismark(reportFiles, colData = data.frame(row.names = ids), strandCollapse = FALSE,
                      nThread=1, BACKEND = "HDF5Array", dir=tempDir, BPPARAM = MulticoreParam(opt$cores), verbose=TRUE)
print('Read.bismark completed.')
gc()

## Get CX context
print('Constructing GRanges from first report and attaching to the BSobj...')
context <- fread(reportFiles[1], colClasses = c('factor', 'numeric',
    'factor', 'integer', 'integer', 'factor', 'factor'))
context_gr <- GRanges(seqnames = chrom, IRanges(start = context[[2]], width = 1), strand = context[[3]],
                      c_context = Rle(context[[6]]), trinucleotide_context = Rle(context[[7]]))
  
## Re-order based on strand
context_gr <- c(context_gr[strand(context_gr) == '+'], context_gr[strand(context_gr) == '-'])
stopifnot(identical(ranges(granges(BSobj)), ranges(context_gr))) 
rowRanges(BSobj) <- context_gr
rm(context, context_gr)
print('Done.')
  
#################################################################################################
#  Split objects by context and save everything to disk
#################################################################################################

#    Write the filtered bsseq objects (the GRanges being the primary attribute written)
print('Splitting BSobj by cytosine context...')
BS_CpG = BSobj[which(rowRanges(BSobj)$c_context == 'CG'),]
BS_nonCpG = BSobj[which(rowRanges(BSobj)$c_context != 'CG'),]

#    Write the methylation and coverage assays (realized into memory in chunks)
print('Realizing the CpG assays and writing to .h5 file...')
M_CpG = writeHDF5Array(assays(BS_CpG)$M, assayPath, 'M_CpG', verbose=TRUE)
Cov_CpG = writeHDF5Array(assays(BS_CpG)$Cov, assayPath, 'Cov_CpG', verbose=TRUE)
print('Realizing the nonCpG assays and writing to .h5 file...')
M_nonCpG = writeHDF5Array(assays(BS_nonCpG)$M, assayPath, 'M_nonCpG', verbose=TRUE)
Cov_nonCpG = writeHDF5Array(assays(BS_nonCpG)$Cov, assayPath, 'Cov_nonCpG', verbose=TRUE)
rm(BSobj)
gc()

#    Point split objects to their new assays on disk
assays(BS_CpG)$M = M_CpG
assays(BS_CpG)$Cov = Cov_CpG
assays(BS_nonCpG)$M = M_nonCpG
assays(BS_nonCpG)$Cov = Cov_nonCpG

#    Save the BSobjs (minus their assays)
print('Writing the remaining components of the BSobjs...')
save(BS_CpG, file=paste0('/bs_', opt$chr, '_CpG.rda'))
save(BS_nonCpG, file=paste0('/bs_', opt$chr, '_nonCpG.rda'))
print('Done with all tasks.')

#    Ensure R doesn't save a .RData file, potentially on the order of 10s of GBs
rm(list = ls())
