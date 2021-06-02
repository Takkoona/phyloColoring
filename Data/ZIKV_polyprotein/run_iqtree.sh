#!/usr/bin/env bash

samplingTimes=1000

iqtree -pre gtr_ufboot$samplingTimes \
       -m LG+I+G4 \
       -s aligned.fasta \
       -o $(cat ./outgroup.txt) \
       -bb $samplingTimes \
       -bnni \
       -nt AUTO
       
