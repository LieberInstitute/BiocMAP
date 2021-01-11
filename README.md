# WGBS Pipeline

![Workflow Overview](https://github.com/LieberInstitute/WGBS-Pipeline/blob/master/workflow.png)

This pipeline consists of two "modules", which together take a set of FASTQ files, described in a `samples.manifest` file, and ultimately produce [bsseq](https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html) R objects containing methylation data and an R data frame of alignment and quality metrics.

We provide tips for achieving significant increases in throughput and customization [LINK TO MANUSCRIPT WHEN IT'S OUT], by implementing the earlier processing steps manually in place of the "first module" we provide in this repository. We recommend this manual approach to advanced users who handle large WGBS datasets or are particularly interested in performance. Otherwise, one can run the first module and second module in series for a complete WGBS processing workflow.

## First Module

### Inputs

This module processes WGBS reads in FASTQ format. Data may be paired-end or single-end, and a single sample may be split across multiple FASTQ files. A `samples.manifest` file describes the samples to be processed by the pipeline (see the `--input` flag). The `samples.manifest` file associates each FASTQ file with a path and ID, and allows the pipeline to automatically merge files if necessary. Each line in `samples.manifest` should have the following format:

+ *For a set of unpaired reads* `<PATH TO FASTQ FILE>(tab)<optional MD5>(tab)<sample label/id>`
+ *For paired-end sets of reads* `<PATH TO FASTQ 1>(tab)<optional MD5 1>(tab)<PATH TO FASTQ 2>(tab)<optional MD5 2>(tab)<sample label/id>`

A line of paired-end reads could look like this:

`WGBS_sample1_read1.fastq    0    WGBS_sample1_read2.fastq    0    sample1`

### Outputs

At the end of the processing done by the first module, deduplicated and quality-filtered alignments in SAM format are produced for each sample. A `rules.txt` file is also produced, a required input to the second processing module.

## Second Module

### Inputs

**A `samples.manifest` file**, described in the section on inputs to the first module, is also a required input for the second module. For those running the first and second module in series, use the same `samples.manifest` file.

**A `rules.txt` file**, describing where relevant logs and inputs are located. After running the first module, `rules.txt` is automatically created in the directory where `samples.manifest` was read.

```
# An example 'rules.txt' file
manifest = /users/nick/samples.manifest
sam = /users/nick/Arioc/[id]/sams/[id].cfu.sam
arioc_log = /users/nick/Arioc/[id]/log/AriocP.[id].log
bme_log = /users/nick/Arioc/[id]/log/[id].cfu.BME.log
xmc_log = /users/nick/Arioc/[id]/log/[id].cfu.XMC.log
trim_report = /users/nick/trim_galore/[id]/[id].fastq.gz_trimming_report.txt
```

- This file consists of several lines of key-value pairs, where keys and values are separated by an equals sign ("=") and optionally spaces.
- Lines without "=" are ignored, and can be used as comments. The above example starts such lines with "#" for clarity.
- The required keys to include in a valid `rules.txt` file include "manifest", "sam", "arioc_log", and "trim_report". The associated values for these keys are the paths to the `samples.manifest` file, the filtered/deduplicated alignment SAMs, verbose output logs from Arioc alignment, and output logs from `TrimGalore!`, respectively.
- Optional keys accepted in a `rules.txt` file include "bme_log" and "xmc_log". Associated values are paths to logs from Bismark Methylation Extractor and logs from XMC [CLARIFICATION NEEDED!], respectively. Typically, methylation extraction is done in the second module, and these options are only included for compatibility with a previous WGBS processing approach. Information from these logs will be included in the R data frame output from the second module, if these logs are present and specified in `rules.txt`.
- Since exact paths are different between samples, including "[id]" in the value field in a line of `rules.txt` indicates that the path for any particular sample can be found by replacing "[id]" with its sample ID. Otherwise, paths are interpreted literally.
- The `--input [dir]` argument to the `nextflow` command in `run_second_half_*.sh` scripts specifies the directory containing the `rules.txt` file, and is a required argument.

### Outputs

The major outputs from the second module are R objects from the [Bioconductor](https://bioconductor.org/) package (https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html), which contain methylation proportion and coverage information at all cytosine loci in the human genome. `bsseq` extends the [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) class, which provides a general and popular format for storing genomics data. Two `bsseq` objects are produced, with one object containing cytosine sites in CpG context, and the other containing the remaining CpH loci.

#### CpG bsseq object

We "strand collapse" CpG loci, which involves combining methylation data from both genomic strands (and thus discarding strand-specific information). The object is "smoothed" with the [BSmooth algorithm](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r83), a process for inferring  stable regional estimates of methylation levels. Loci are ordered by genomic position.

#### CpH bsseq object

We retain strand-specific information for CpH loci. These loci are also ordered by genomic position.

#### Storage method

Because all cytosines in the genome are included in the output objects, the data may occupy tens or even hundreds of gigabytes in memory (RAM) if loaded in a typical fashion. To enable working with the objects in a reasonable amount of memory, the assays (in this case methylation fraction and coverage counts) are HDF5-backed using the [HDF5Array](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html) R package. Essentially, this involves storing the assays in a `.h5` file, a format designed to enable working with on-disk data as if it were loaded in RAM.

#### R data frame

[TODO!!] rows are sample names, and columns must be individually described here.
