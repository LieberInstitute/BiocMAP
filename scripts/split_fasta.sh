#!/bin/bash

FASTA=$1

numChrs=$(grep ">c" $FASTA | wc -l)
lineNums=$(grep -n ">" $FASTA | cut -d : -f 1 | paste -sd " ")
for CHR in `seq 1 $numChrs`; do
    start=$(echo $lineNums | cut -d " " -f $CHR)
    end=$(echo $lineNums | cut -d " " -f $(($CHR+1)))
    chr_name=$(sed -n "${start}p;$(($start+1))q" $FASTA | cut -d ">" -f 2 | cut -d " " -f 1)
    sed -n "$start,$(($end-1))p;$endq" $FASTA > ${chr_name}.fa
done