#  This script exists to document a semi-automated method for producing a list
#  of R packages, read by "checkpoint::checkpoint" when installing R and
#  required packages via the "local" or "conda" installation modes. At the time
#  of producing this script, "checkpoint::checkpoint" halts with an error if
#  any Bioconductor packages are present in the directory passed to
#  "project_dir", thus motivating a workaround where a dummy directory
#  including only non-Bioc library calls is used.

grep -E "library\(.*\)" ../*.R \
    | cut -d ":" -f 2 \
    | tr  -d " " \
    | sed 's/suppressPackageStartupMessages(\(.*\))/\1/' \
    | sort -u > r_package_list.R
    
#  Then manually remove lines associated with Bioconductor packages and
#  "jaffelab", and add "BiocManager" and "remotes"