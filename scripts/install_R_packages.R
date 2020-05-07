ordinary_packages = c('devtools', 'getopt', 'data.table', 'BiocManager')
bioc_packages = c('bsseq', 'GenomicRanges', 'HDF5Array')

install.packages(ordinary_packages)

library('BiocManager')
BiocManager::install(bioc_packages)

library('devtools')
install_github('LieberInstitute/jaffelab')
