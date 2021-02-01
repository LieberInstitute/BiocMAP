#!/bin/bash

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/bin/nextflow first_half.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile first_half_local \
    > run_first_half_local.log
