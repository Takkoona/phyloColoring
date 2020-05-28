#!/usr/bin/python3
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

PDB_DIR = "."
MODELLING_DIR = os.path.join(PDB_DIR, "Modelling")

TEMPLATE_ID = "4lxv"
PDB_PARSER = PDBParser()

def calculateRMS(row):
    rmsRow = []
    for query, subject in row:
        if query.id == subject.id:
            rmsRow.append(0.0)
        else:
            x = np.array([res["CA"].coord for res in query.get_residues()])
            y = np.array([res["CA"].coord for res in subject.get_residues()])
            sup = SVDSuperimposer()
            sup.set(x, y)
            sup.run()
            rmsRow.append(sup.get_rms())
    return rmsRow

if __name__ == '__main__':
    structures = []
    # n = 0

    for seqDir in os.listdir(MODELLING_DIR):
        seqName = os.path.basename(seqDir)
        structure = None
        align = None
        for pdbFile in glob(os.path.join(MODELLING_DIR, seqDir, "*.pdb")):
            structure = PDB_PARSER.get_structure(seqName, pdbFile)
        for alignFile in glob(os.path.join(MODELLING_DIR, seqDir, seqName + "-" + TEMPLATE_ID +".ali")):
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
        s = PDB_PARSER.get_structure(seqName, f)
        structures.append(s)
            
        # n += 1
        # if n == 8:
        #     break

    print("Start rms calculation...")

    p = Pool(cpu_count())
    res = p.map(
        calculateRMS, 
        ([(query, subject) for subject in structures] for query in structures)
    )
    seqNames = [s.id for s in structures]
    res = pd.DataFrame(res, columns=seqNames, index=seqNames)
    res.to_csv(os.path.join(PDB_DIR, "rms_matrix.csv"))
