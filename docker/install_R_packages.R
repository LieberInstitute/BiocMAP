#  While the main version of this script is in the 'scripts' directory, a
#  specific copy used while producing the 'bioc_kallisto:3.11' docker image is
#  kept here for reproducibility, in case changes are made to the original
#  script.

ordinary_packages = c('remotes', 'getopt', 'data.table', 'BiocManager', 'here')
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
