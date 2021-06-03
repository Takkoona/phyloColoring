#!/usr/bin/env bash

iqtree -pre gtr_ufboot \
       -m GTR20+I+G4 \
       -s HA1_aligned.fasta \
       -o AF201874.1 \
       -bb 1000 \
       -bnni \
       -nt AUTO
