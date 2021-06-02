#!/usr/bin/env bash

iqtree -pre gtr_ufboot \
       -m GTR20+I+G4 \
       -s aligned.fasta \
       -o HK_1_1968 \
       -bb 1000 \
       -bnni \
       -nt AUTO
