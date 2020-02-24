library('data.table')
library('GenomicRanges')
library('getopt')

spec <- matrix(c(
    'threads', 't', 1, 'integer', 'number of threads to utilize'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

data.table::setDTthreads(opt$threads)

#  Read cytosine report into memory
context = fread(list.files(pattern='*.cytosine_report.txt'), 
                colClasses = c('factor', 'numeric','factor', 'integer', 
                               'integer', 'factor', 'factor'))
       
for (c_context in c("CG", "CHG", "CHH")) {
    temp = which(context[,6] == c_context)
    if (length(temp) > 0) {
        m = sum(context[temp, 4])
        u = sum(context[temp, 5])
        print(paste0("C methylated in ", c_context, " context: ", round(100 * m/ (m + u), 1), "%"))
    } else {
        print(paste0("C methylated in ", c_context, " context: undefined (no sites)"))
    }
}
