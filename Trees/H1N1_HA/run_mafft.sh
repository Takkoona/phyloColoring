#!/bin/bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N H1N1_HA_MSA

input="sequences.fasta"
output="aligned.fasta"

mafft --auto $input > $output
