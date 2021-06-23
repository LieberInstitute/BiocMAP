ordinary_packages = c('remotes', 'getopt', 'data.table', 'BiocManager')
bioc_packages = c('bsseq', 'GenomicRanges', 'HDF5Array', 'BiocParallel')

#  Ordinary R packages
for (p in ordinary_packages) {
    if (!requireNamespace(p, quietly=TRUE)) {
        install.packages(p)
    }
}

#  Bioconductor packages
for (p in bioc_packages) {
    if (!requireNamespace(p, quietly=TRUE)) {
        BiocManager::install(p)
    }
}

#  GitHub packages
if (!requireNamespace('jaffelab', quietly=TRUE)) {
    remotes::install_github('LieberInstitute/jaffelab')
}
