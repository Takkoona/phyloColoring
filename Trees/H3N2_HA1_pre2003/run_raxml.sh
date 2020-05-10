#!/bin/bash

input="aligned.fasta"
outgroup="HK/1/1968"

rm *.PROTGAMMAGTR

raxmlHPC -s $input \
         -o $outgroup \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -p 12345
