#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N ZIKV_MSA

input="aligned.fasta"
outgroup="ANK57896"

rm *.PROTGAMMAGTR

raxmlHPC -s $input \
         -o $outgroup \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -p 12345
