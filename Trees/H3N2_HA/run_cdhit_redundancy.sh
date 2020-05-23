#!/bin/bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H3N2_HA_redundancy

dirPath="Spatiotemporal"

for input in $dirPath/*/outliered.fasta
do
    output=${input//outliered/trimmed}
    cdhit -i "$input" -o "$output" -c 0.998
done
