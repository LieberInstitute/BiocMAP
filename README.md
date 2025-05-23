[![DOI](https://zenodo.org/badge/223444482.svg)](https://zenodo.org/badge/latestdoi/223444482)

# BiocMAP

BiocMAP is a **Bioc**onductor-friendly **M**ethylation **A**nalysis **P**ipeline. It consists of two [nextflow](https://www.nextflow.io/)-based "modules", which together take a set of FASTQ files, described in [a `samples.manifest` file](http://research.libd.org/BiocMAP/inputs.html#the-samples.manifest-file), and ultimately produce [bsseq](https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html) R objects containing methylation data and an R data frame of alignment and quality metrics.

The first BiocMAP module performs speedy alignment to a reference genome by [Arioc](https://github.com/RWilton/Arioc), and requires GPU resources. Methylation extraction and remaining steps are performed in the second module, optionally on a different computing system where GPUs need not be available.

![Workflow Overview](https://github.com/LieberInstitute/BiocMAP/blob/master/workflow.png)

## Features

- GPU-accelerated alignment to a reference genome via [Arioc](https://github.com/RWilton/Arioc)
- Memory-efficient, HDF5-backed [bsseq](https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html) output objects immediately ready for analysis with [Bioconductor](https://bioconductor.org/)/R packages of choice
- Automatic management of [reference files](http://research.libd.org/BiocMAP/annotation.html), allowing simple configuration of [GENCODE](https://www.gencodegenes.org/) release and [other settings](http://research.libd.org/BiocMAP/annotation.html#choosing-build), while alternatively [supporting user-provided files](http://research.libd.org/BiocMAP/annotation.html#custom-annotation)
- [Support for docker and singularity](http://research.libd.org/BiocMAP/setup-details.html#installation) for flexible and reproducible installation
- Automatically merge samples split across multiple FASTQ files, using [the `samples.manifest` input](http://research.libd.org/BiocMAP/inputs.html#the-samples.manifest-file)

## Get started

[The BiocMAP documentation website](http://research.libd.org/BiocMAP/index.html) provides a complete description of features, installation, and many other details. To quickly get started, see our [quick-start guide](http://research.libd.org/BiocMAP/quick-start.html).

We provide [shell scripts](http://research.libd.org/BiocMAP/quick-start.html#your-main-script) for out-of-the-box execution a [SLURM](https://slurm.schedmd.com/overview.html) or [SGE](https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html)-based computing cluster, or for execution on a Linux-based machine. Software dependencies can be installed via the shell script `install_software.sh`, which makes use of docker or Anaconda/Miniconda.

## Authors

[Nick Eagles](https://github.com/Nick-Eagles)

## Cite

We hope `BiocMAP` will be a useful tool for your research. Please use the following bibtex information to cite this software. Thank you!

```
@article {Eagles2023,
	author = {Eagles, Nicholas J and Wilton, Richard and Jaffe, Andrew E. and Collado-Torres, Leonardo},
	title = {BiocMAP: A Bioconductor-friendly, GPU-accelerated pipeline for bisulfite-sequencing data},
	year = {2023},
	doi = {10.1186/s12859-023-05461-3},
	publisher = {Springer Science and Business Media LLC},
	URL = {https://doi.org/10.1186/s12859-023-05461-3},
	journal = {BMC Bioinformatics}
}
```

## Contact

[Leonardo Collado-Torres](http://lcolladotor.github.io/)
