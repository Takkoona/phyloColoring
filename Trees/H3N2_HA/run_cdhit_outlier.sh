#!/bin/bash

dirPath="Spatiotemporal"

for input in $dirPath/*/sequences.fasta
do
    output=${input//sequences/clustered}
    cdhit -i "$input" -o "$output" -c 0.95
done
