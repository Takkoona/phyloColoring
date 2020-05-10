#!/bin/bash

dirPath="Spatiotemporal"

for input in $dirPath/*/outliered.fasta
do
    output=${input//outliered/trimmed}
    cdhit -i "$input" -o "$output" -c 0.998
done
