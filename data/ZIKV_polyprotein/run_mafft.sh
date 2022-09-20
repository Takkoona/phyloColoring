#!/usr/bin/env bash

input="sequences.fasta"
output="aligned.fasta"

mafft --auto $input > $output
