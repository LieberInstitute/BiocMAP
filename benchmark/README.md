This directory contains all code used to benchmark BiocMAP against comparable alternatives, as part of the manuscript.

The following directories include code to install and run different pipeline tools on the example data:

Published as part of the manuscript benchmark:

- **BiocMAP_dir**: [BiocMAP](https://github.com/LieberInstitute/BiocMAP)
- **methylseq**: [nf-core/methylseq](https://github.com/nf-core/methylseq)

Tested, but did not make the manuscript due to errors or prohibitively long execution time:

- **bs_seeker3**: [BS-Seeker3](https://github.com/khuang28jhu/bs3)
- **wg_blimp**: [wg-blimp](https://github.com/MarWoes/wg-blimp)

In the `benchmark_stats` directory, custom R code and the `sgejobs` package were used to trace the jobs run for each pipeline, as reported by the SGE job scheduler, and compile information about run time and computational resource usage.
