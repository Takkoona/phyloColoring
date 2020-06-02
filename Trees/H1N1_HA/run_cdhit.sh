#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H1N1_HA_trim

dirPath="Spatiotemporal"

for input in $dirPath/*/sequences.fasta
do
    output=${input//sequences/trimmed}
    cdhit -i "$input" -o "$output" -c 0.998
done
