#  While the main version of this script is in the 'scripts' directory, a
#  specific copy used while producing the 'bioc_kallisto:3.17' docker image is
#  kept here for reproducibility. In addition, packages not strictly required
#  by this pipeline are installed into the image so that the image can also be
#  used by SPEAQeasy (https://github.com/LieberInstitute/SPEAQeasy).

ordinary_packages = c(
    'remotes', 'getopt', 'data.table', 'BiocManager', 'here', 'devtools',
    'matrixStats', 'plyr', 'rafalib', 'RColorBrewer', 'sessioninfo', 'styler'
)

bioc_packages = c(
    'bsseq', 'GenomicRanges', 'HDF5Array', 'BiocParallel', 'biomaRt',
    'Biostrings', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
    'BSgenome.Mmusculus.UCSC.mm10', 'BSgenome.Rnorvegicus.UCSC.rn6',
    'DelayedArray', 'derfinder', 'GenomicFeatures', 'org.Hs.eg.db',
    'org.Mm.eg.db', 'org.Rn.eg.db', 'rtracklayer', 'SummarizedExperiment'
)

#  Ordinary R packages
for (p in ordinary_packages) {
    if (!requireNamespace(p, quietly=TRUE)) {
        install.packages(p, repos='http://cran.us.r-project.org')
    }
}

#  Bioconductor packages
for (p in bioc_packages) {
    if (!requireNamespace(p, quietly=TRUE)) {
        BiocManager::install(p, update = FALSE)
    }
}

#  GitHub packages: the version of GenomeInfoDb that comes with Bioc 3.17
#  has this issue: https://github.com/Bioconductor/GenomeInfoDb/issues/86,
#  motivating the explicit reinstall of a fixed version here
remotes::install_github('LieberInstitute/jaffelab')
remotes::install_github('Bioconductor/GenomeInfoDb@284ef2a')
