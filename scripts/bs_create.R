#  Read cytosine reports and produce 2 HDF5-backed BSseq objects (CpG and nonCpG context).

suppressPackageStartupMessages(library('bsseq'))
suppressPackageStartupMessages(library('HDF5Array'))
suppressPackageStartupMessages(library('BiocParallel'))
suppressPackageStartupMessages(library('getopt'))
suppressPackageStartupMessages(library('GenomicRanges'))
suppressPackageStartupMessages(library('jaffelab'))
library('data.table')

r_time = proc.time()

################################################################
#  Convenience variables
################################################################

spec <- matrix(c(
    'chr', 's', 1, 'character', 'chromosome to subset',
    'region', 'r', 1, 'character', 'brain region',
    'cores', 'c', 1, 'integer', 'number of cores to utilize',
    'dir', 'd', 1, 'character', 'output dir for objects'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

tempDir = paste0('temp_hd5/', opt$chr, '/')
dir.create('temp_hd5')

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

reportFiles = list.files(pattern='.*\\.CX_report\\.txt')
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
context_gr <- GRanges(seqnames = opt$chr, IRanges(start = context[[2]], width = 1), strand = context[[3]],
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

#  Note the entire CpG object is realized into memory- at the time of writing
#  this script, strand collapsing takes an unreasonable time to block process
print('Subsetting to CpG context, strand collapsing, and saving...')
dir.create(
    file.path(opt$dir, opt$chr, 'CpG'), recursive=TRUE, showWarnings = FALSE
)
BS_CpG = strandCollapse(realize(BSobj[which(rowRanges(BSobj)$c_context == 'CG'),]))
BS_CpG = saveHDF5SummarizedExperiment(BS_CpG,
                                      dir=file.path(opt$dir, opt$chr, 'CpG'),
                                      replace=TRUE)
gc()

#  Perform smoothing and assign to dummy object
print('Smoothing HDF5-backed CpG object...')
bs_smooth = BSmooth(BS_CpG, 
                    BPPARAM = MulticoreParam(opt$cores), 
                    verbose = TRUE)

rm(BS_CpG, bs_smooth)
gc()

print('Subsetting to CpH context, sorting ranges, and saving...')
dir.create(
    file.path(opt$dir, opt$chr, 'CpH'), recursive=TRUE, showWarnings = FALSE
)
BS_CpH = BSobj[which(rowRanges(BSobj)$c_context != 'CG'),]
BS_CpH = BS_CpH[order(ranges(BS_CpH)),]
saveHDF5SummarizedExperiment(BS_CpH, dir=file.path(opt$dir, opt$chr, 'CpH'),
                             replace=TRUE)

print('Done with all tasks.')

gc()
proc.time() - r_time

#    Ensure R doesn't save a .RData file, potentially on the order of 10s of GBs
rm(list = ls())
