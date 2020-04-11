#!/bin/bash

input="sequences.fasta"
output="aligned.fasta"

mafft $input > $output
