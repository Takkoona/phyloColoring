#!/bin/bash

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
