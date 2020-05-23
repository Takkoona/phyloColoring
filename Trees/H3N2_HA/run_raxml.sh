#!/bin/bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H3N2_HA_tree

input="aligned.fasta"
outgroup="AB284320"

if ls . | grep .PROTGAMMAGTR
then
    rm *.PROTGAMMAGTR
fi

raxmlHPC -s $input \
         -o $outgroup \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -p 12345
