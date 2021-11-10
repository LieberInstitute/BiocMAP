# BiocMAP

BiocMAP is a **Bioc**onductor-friendly **M**ethylation **A**nalysis **P**ipeline. It consists of two [nextflow](https://www.nextflow.io/)-based "modules", which together take a set of FASTQ files, described in [a `samples.manifest` file](http://research.libd.org/WGBS-Pipeline/inputs.html#the-samples.manifest-file), and ultimately produce [bsseq](https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html) R objects containing methylation data and an R data frame of alignment and quality metrics.

We provide tips for achieving significant increases in throughput and customization [LINK TO MANUSCRIPT WHEN IT'S OUT], by implementing the earlier processing steps manually in place of the "first module" we provide in this repository. We recommend this manual approach to advanced users who handle large WGBS datasets or are particularly interested in performance. Otherwise, one can run the first module and second module in series for a complete WGBS processing workflow.

![Workflow Overview](https://github.com/LieberInstitute/WGBS-Pipeline/blob/master/workflow.png)

## Features

- GPU-accelerated alignment to a reference genome via [Arioc](https://github.com/RWilton/Arioc)
- Memory-efficient, HDF5-backed [bsseq](https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html) output objects immediately ready for analysis with [Bioconductor](https://bioconductor.org/)/R packages of choice
- Automatic management of [reference files](http://research.libd.org/WGBS-Pipeline/annotation.html), allowing simple configuration of [GENCODE](https://www.gencodegenes.org/) release and [other settings](http://research.libd.org/WGBS-Pipeline/annotation.html#choosing-build), while alternatively [supporting user-provided files](http://research.libd.org/WGBS-Pipeline/annotation.html#custom-annotation)
- [Support for docker and conda](http://research.libd.org/WGBS-Pipeline/setup-details.html#installation) for flexible and reproducible installation
- Automatically merge samples split across multiple FASTQ files, using [the `samples.manifest` input](http://research.libd.org/WGBS-Pipeline/inputs.html#the-samples.manifest-file)

## Get started

[The BiocMAP documentation website](http://research.libd.org/WGBS-Pipeline/index.html) provides a complete description of features, installation, and many other details. To quickly get started, see our [quick-start guide](http://research.libd.org/WGBS-Pipeline/quick-start.html).

We provide [shell scripts](http://research.libd.org/WGBS-Pipeline/quick-start.html#your-main-script) for out-of-the-box execution a [SLURM](https://slurm.schedmd.com/overview.html) or [SGE](https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html)-based computing cluster, or for execution on a Linux-based machine. Software dependencies can be installed via the shell script `install_software.sh`, which makes use of docker or Anaconda/Miniconda.

## Authors

[Nick Eagles](https://github.com/Nick-Eagles)

## Contact

[Leonardo Collado-Torres](http://lcolladotor.github.io/)
