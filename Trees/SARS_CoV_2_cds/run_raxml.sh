#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N SARS-CoV-2_MSA

input="aligned.fasta"
#outgroup="$(cat outgroup.txt)"
outgroup="EPI_ISL_402125"

rm *.GTRGAMMA

raxmlHPC -s $input \
         -o $outgroup \
         -n GTRGAMMA \
         -m GTRGAMMA \
         -p 12345
