#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q NGS
#$ -pe mpi 24
#$ -cwd
#$ -N SARS-CoV-2_MSA

input="sequences.fasta"
output="aligned.fasta"

mafft $input > $output
