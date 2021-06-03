#!/usr/bin/env bash

raxmlHPC -f a \
         -p 12345 \
         -x 12345 \
         -# 1000 \
         -s HA1_aligned.fasta \
         -n PROTGAMMAGTR \
         -m PROTGAMMAGTR \
         -o AF201874.1
