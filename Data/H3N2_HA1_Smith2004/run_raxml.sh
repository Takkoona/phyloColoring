#!/usr/bin/env bash

raxmlHPC -f a \
         -p 12345 \
         -x 12345 \
         -# 1000 \
         -s aligned.fasta \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -o HK/1/1968
