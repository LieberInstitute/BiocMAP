#  This script is only intended to be run from 'install_software.sh', where the
#  working directory is [BiocMAP dir]/Software

print("Checking R packages...")

for (package in c("checkpoint", "here")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, repos = "http://cran.us.r-project.org")
    }
}

library("here")
library("checkpoint")

#  Install ordinary packages as they existed when R 4.1.0 was released
dir.create(here('Software', '.checkpoint'))
checkpoint("2021-09-01",
    project_dir = here("scripts", "r_packages"),
    checkpoint_location = here('Software')
)

#  Install the Bioc packages and the GitHub package "jaffelab"
bioc_packages = c(
    'bsseq', 'GenomicRanges', 'HDF5Array', 'BiocParallel',
    'LieberInstitute/jaffelab'
)

BiocManager::install(packages, update = FALSE)
