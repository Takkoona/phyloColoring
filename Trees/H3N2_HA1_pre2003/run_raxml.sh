#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N Smith2004_tree

input="aligned.fasta"
outgroup="HK/1/1968"

rm *.PROTGAMMAGTR

raxmlHPC -s $input \
         -o $outgroup \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -p 12345
