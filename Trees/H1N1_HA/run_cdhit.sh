#!/bin/bash

dirPath="Spatiotemporal"

for input in $dirPath/*/sequences.fasta
do
    output=${input//sequences/trimmed}
    cdhit -i "$input" -o "$output" -c 0.998
done
