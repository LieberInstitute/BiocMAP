---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
title: "[BUG] Your bug or feature request"
labels: ''
assignees: ''
---

Please briefly describe your problem and what output/behavior you expect.

## Context

Provide some context for your bug report or feature request. This could be the:

* link to raw code, example: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/00-template.Rmd#L24-L28
* link to a commit, example: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/commit/6aa30b22eda614d932c12997ba611ba582c435d7
* link to a line of code inside a commit, example: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/commit/6aa30b22eda614d932c12997ba611ba582c435d7#diff-e265269fe4f17929940e81341b92b116R17
* link to code from an R package, example: https://github.com/LieberInstitute/spatialLIBD/blob/master/R/run_app.R#L51-L55

## Code

Include the code you ran and comments

```R
## prompt an error
stop('hola')

## check the error trace
traceback()
```

## Error message

When reporting a bug or error, please include the full error message reported in the main output log (e.g. for SLURM users and the first module, `run_first_half_slurm.log`). This may be short, but should (usually) include the process name where the error occurred, and a follow-up explanation:

```
Error executing process > 'EncodeReference'

Caused by:
  Cannot get property 'annotation' on null object
```


## BiocMAP run information

For errors or issues during BiocMAP runs, please copy and paste the header printed in the main output log (e.g. for SLURM users and the first module, `run_first_half_slurm.log`). The header will look similar to this:

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

