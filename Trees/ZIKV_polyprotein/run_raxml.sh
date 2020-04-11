#!/bin/bash

input="aligned.fasta"
outgroup="ANK57896"

rm *.PROTGAMMAGTR

raxmlHPC -s $input \
         -o $outgroup \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -p 12345
