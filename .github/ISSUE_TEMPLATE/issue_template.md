---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
title: "[BUG] Your bug or feature request"
labels: ''
assignees: ''
---

## Feature requests

For feature requests, please describe the desired functionality and how it differs from what's currently implemented. If the benefits for the suggested feature aren't obvious, please provide a brief explanation motivating the change.

The following example is based off of [issue 14](https://github.com/LieberInstitute/BiocMAP/issues/14):

### Current functionality

The `Trimming` subdirectory in the output directory contains uncompressed FASTQ files, even for samples that aren't trimmed.

### Desired functionality

Only samples that undergo trimming should produce trimmed FASTQ files that are published to the output directory. Untrimmed FASTQ files should never be published as outputs.

### Motivation for change

Untrimmed FASTQs are duplicate data, and are confusing to include in the "Trimming" directory. The FASTQs are also uncompressed, meaning a great portion of disk space is unnecessarily taken for these temporary files.

## Bugs/ error reports

### Error message

When reporting a bug or error, please include the full error message reported in the main output log (e.g. for SLURM users and the first module, `run_first_half_slurm.log`). This may be short, but should (usually) include the process name where the error occurred, and a follow-up explanation:

```
Error executing process > 'EncodeReference'

Caused by:
  Cannot get property 'annotation' on null object
```


### BiocMAP run information

Please copy and paste the header printed in the main output log (e.g. for SLURM users and the first module, `run_first_half_slurm.log`). The header will look similar to this:

```
================================================================================
    BiocMAP- First Module
================================================================================
---- Main options:
BiocMAP version    = 60766c92528f1f302fca1350179131ad0742a741
Config profile     = first_half_slurm
All alignments     = false
Annotation dir     = /path/to/annotation
Annotation release = 34
Annotation build   = main
Custom anno label  = 
Input dir          = /input/dir
Output dir         = /output/dir
Reference          = hg38
Sample	           = paired
Trim mode          = force
Working dir        = /work/dir
Current user	   = nickeagles
---- Software arguments:
Arioc GPU batch size  = 32k
Arioc gapped seed     = hsi25_0_30_CT
Arioc non-gapped seed = ssi84_2_30_CT
Arioc gapped args     = Wmxgs="2,6,5,3" Vt="L,0,1" maxJ="20" seedDepth="4"
Arioc non-gapped args = maxJ="200" maxMismatches="5"
Arioc X args          = watchdogInterval="60" cgaReserved="24M" useHinGmem="0" useJinGmem="0" useHJinGPmem="0" serialLUTinit="1"
Arioc max GPUs        = 1
Manually set GPU      = false
GPU usage cutoff      = 10
================================================================================
```

Other helpful details include the installation method originally used to set up BiocMAP (e.g. `singularity`) and the approximate number of samples in your dataset.
