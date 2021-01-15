ordinary_packages = c('devtools', 'getopt', 'data.table', 'BiocManager')
bioc_packages = c('bsseq', 'GenomicRanges', 'HDF5Array', 'BiocParallel')

#  Ordinary R packages
for (p in ordinary_packages) {
    if (!requireNamespace(p, quietly=TRUE)) {
        install.packages(p)
    }
}

#  Bioconductor packages
library('BiocManager')
BiocManager::install(bioc_packages)

#  GitHub packages
library('devtools')
install_github('LieberInstitute/jaffelab')
