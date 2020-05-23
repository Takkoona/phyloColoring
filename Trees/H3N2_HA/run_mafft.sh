#!/bin/bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H3N2_HA_MSA

input="sequences.fasta"
output="aligned.fasta"

mafft --auto $input > $output
