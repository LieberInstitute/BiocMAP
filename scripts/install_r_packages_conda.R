#  This script is only intended to be run from 'install_software.sh "conda"'.
#  Here, we install non-Bioc R packages (with the exception of 'here',
#  'checkpoint',  and 'jaffelab') using MRAN snapshots (the 'checkpoint'
#  package). This ensures packages install properly (and the same way) even
#  when R 4.1.0 becomes "outdated". During development, this seemed less buggy
#  than the theoretically preferable method of installing ordinary R packages
#  through the 'r' conda channel using fixed versions. However, installing
#  Bioconductor packages via the 'bioconda' conda channel seemed to work fine
#  (and 'checkpoint' can't handle Bioconductor packages anyway), so this is
#  performed directly in the 'install_software.sh' script.

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

BiocManager::install('LieberInstitute/jaffelab', update = FALSE)
