#!/bin/bash

FASTA=$1

numChrs=$(grep ">c" $FASTA | wc -l)
lineNums=$(grep -n ">" $FASTA | cut -d : -f 1 | paste -sd " ")
for CHR in $(seq 1 $numChrs); do
    start=$(echo $lineNums | cut -d " " -f $CHR)
    chr_name=$(sed -n "${start}p;$(($start+1))q" $FASTA | cut -d ">" -f 2 | cut -d " " -f 1)
    if [ $CHR == $numChrs ]; then
        #  Range is from start to EOF
        sed -n "$start,$ p;$endq" $FASTA > ${chr_name}.fa
    else
        end=$(echo $lineNums | cut -d " " -f $(($CHR+1)))
        sed -n "$start,$(($end-1))p;$endq" $FASTA > ${chr_name}.fa
    fi
done
