#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H1N1_HA_tree

input="aligned.fasta"
outgroup="FJ966082"

if ls . | grep .PROTGAMMAGTR
then
    rm *.PROTGAMMAGTR
fi

raxmlHPC -s $input \
         -o $outgroup \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -p 12345
