#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -S /usr/bin/python3
#$ -q NGS
#$ -pe mpi 1
#$ -N modeller_alignment

import os
import json
from glob import glob

import modeller
from Bio import SeqIO

PDB_DIR = "."
MODELLING_DIR = os.path.join(PDB_DIR, "Modelling")

TEMPLATE_ID = "4lxv"

with open(os.path.join(PDB_DIR, "tipNames.json")) as f:
    tipNames = json.load(f)
    
env = modeller.environ()

aln = modeller.alignment(env)

for seqDir in os.listdir(MODELLING_DIR):
    seqName = os.path.basename(seqDir)
    if seqName in tipNames:
        for pdbFile in glob(os.path.join(MODELLING_DIR, seqDir, "*.pdb")):
            m = modeller.model(env, file=pdbFile)
            aln.append_model(m, atom_files=pdbFile, align_codes=tipNames[seqName])

aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file=os.path.join(PDB_DIR, 'distance.mat'))
