#!/bin/bash

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/bin/nextflow second_half.nf \
    --sample "paired" \
    --reference "hg38" \
    -profile second_half_local \
    > run_second_half_local.log
