# WGBS Pipeline

## master branch

![Workflow Overview](https://github.com/LieberInstitute/WGBS-Pipeline/blob/master/workflow.png)

The master branch of the pipeline takes a set of FASTQ files, described in a `samples.manifest` file, and ultimately produces [bsseq](https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html) R objects containing methylation data and an R data frame of alignment and quality metrics.

### Inputs

This pipeline processes WGBS reads in FASTQ format. Data may be paired-end or single-end, and a single sample may be split across multiple FASTQ files. A `samples.manifest` file describes the samples to be processed by the pipeline (see the `--input` flag). The `samples.manifest` file associates each FASTQ file with a path and ID, and allows the pipeline to automatically merge files if necessary. Each line in `samples.manifest` should have the following format:

+ *For a set of unpaired reads* `<PATH TO FASTQ FILE>(tab)<optional MD5>(tab)<sample label/id>`
+ *For paired-end sets of reads* `<PATH TO FASTQ 1>(tab)<optional MD5 1>(tab)<PATH TO FASTQ 2>(tab)<optional MD5 2>(tab)<sample label/id>`

A line of paired-end reads could look like this:

`WGBS_sample1_read1.fastq    0    WGBS_sample1_read2.fastq    0    sample1`

### Outputs

The major outputs from this pipeline are R objects from the [Bioconductor](https://bioconductor.org/) package (https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq.html), which contain methylation proportion and coverage information at all cytosine loci in the human genome. `bsseq` extends the [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) class, which provides a general and popular format for storing genomics data. Two `bsseq` objects are produced, with one object containing cytosine sites in CpG context, and the other containing the remaining CpH loci.

#### CpG bsseq object

We "strand collapse" CpG loci, which involves combining methylation data from both genomic strands (and thus discarding strand-specific information). The object is "smoothed" with the [BSmooth algorithm](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r83), a process for inferring  stable regional estimates of methylation levels. Loci are ordered by genomic position.

#### CpH bsseq object

We retain strand-specific information for CpH loci. These loci are also ordered by genomic position.

#### Storage method

Because all cytosines in the genome are included in the output objects, the data may occupy tens or even hundreds of gigabytes in memory (RAM) if loaded in a typical fashion. To enable working with the objects in a reasonable amount of memory, the assays (in this case methylation fraction and coverage counts) are HDF5-backed using the [HDF5Array](https://bioconductor.org/packages/release/bioc/html/HDF5Array.html) R package. Essentially, this involves storing the assays in a `.h5` file, a format designed to enable working with on-disk data as if it were loaded in RAM.
