#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -S /usr/bin/python3
#$ -q NGS
#$ -pe mpi 24
#$ -N superimpose_mpi

import os
import json
from glob import glob
from io import StringIO
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Dice

MODELLING_DIR = os.path.join(".", "Modelling")
TEMPLATE_ID = "4lxv"

def trimStructure(seqName):
    parser = PDBParser()
    structure = None
    align = None
    for pdbFile in glob(os.path.join(MODELLING_DIR, seqName, "*.pdb")):
        structure = parser.get_structure(seqName, pdbFile)
    for alignFile in glob(os.path.join(MODELLING_DIR, seqName, seqName + "-" + TEMPLATE_ID +".ali")):
        align = SeqIO.index(alignFile, "pir")[TEMPLATE_ID]
        for (start, res) in enumerate(align, start=1):
            if res != '-':
                break
        for (end, res) in enumerate(align[::-1]):
            if res != '-':
                break
    f = StringIO()
    Dice.extract(structure, chain_id=' ', start=start, end=len(align) - end, filename=f)
    f = StringIO(f.getvalue())
    s = parser.get_structure(seqName, f)
    return(s)

def calculateRMS(row):
    scores = {}
    for query, subject in row:
        if query.id == subject.id:
            scores[subject.id] = 0.0
        else:
            x = np.array([res["CA"].coord for res in query.get_residues()])
            y = np.array([res["CA"].coord for res in subject.get_residues()])
            sup = SVDSuperimposer()
            sup.set(x, y)
            sup.run()
            scores[subject.id] = sup.get_rms()
    with open(os.path.join(MODELLING_DIR, query.id, query.id + "RMS.json"), 'w') as f:
        json.dump(scores, f)

if __name__ == '__main__':
    p = Pool(cpu_count())

    structures = p.map(
        trimStructure,
        os.listdir(MODELLING_DIR)[0:6]
    )

    print("Start rms calculation...")

    p.map(
        calculateRMS, 
        ([(query, subject) for subject in structures] for query in structures)
    )
