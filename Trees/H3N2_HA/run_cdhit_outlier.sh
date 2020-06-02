#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H3N2_HA_outlier

dirPath="Spatiotemporal"

for input in $dirPath/*/sequences.fasta
do
    output=${input//sequences/clustered}
    cdhit -i "$input" -o "$output" -c 0.95
done
