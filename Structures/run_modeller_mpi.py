#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -S /usr/bin/python
#$ -q NGS
#$ -pe mpi 24
#$ -N modeller_mpi

import os
from glob import glob
from multiprocessing import Pool, cpu_count

import modeller
from modeller.automodel import automodel, assess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

PDB_DIR = "."
MODELLING_DIR = os.path.join(PDB_DIR, "Modelling")
CURRENT_DIR = os.getcwd()
TEMPLATE_ID = "4lxv"

if not os.path.exists(MODELLING_DIR):
    os.mkdir(MODELLING_DIR)

allSeqs = [record for record in SeqIO.parse(os.path.join(PDB_DIR, "H1N1_HA.fasta"), "fasta")]

def run_modeller(record):
    env = modeller.environ()

    record.description = ':'.join(["sequence", record.id, "", "", "", "", "", "", "0.00", "0.00"])
    record.name = ""
    record.seq.alphabet = generic_protein
    
    outDir = os.path.join(MODELLING_DIR, record.id)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    os.chdir(outDir)
    if glob("*.pdb"):
        os.chdir(CURRENT_DIR)
        return "EXIST"
    
    pdbFile = os.path.join("..", "..", TEMPLATE_ID + ".pdb")
    mdl = modeller.model(env, file=pdbFile, model_segment=('FIRST:A','LAST:B'))
    
    seqFile = record.id + ".ali"
    SeqIO.write(record, seqFile, "pir")
    
    aln = modeller.alignment(env)
    aln.append_model(mdl, align_codes=TEMPLATE_ID, atom_files=pdbFile)
    aln.append(file=seqFile, align_codes=record.id)
    
    aln.align2d()
    alnFile = record.id + "-" + TEMPLATE_ID + ".ali"
    aln.write(file=alnFile, alignment_format='PIR')
    
    a = automodel(
        env=env,
        alnfile=alnFile,
        knowns=TEMPLATE_ID, 
        sequence=record.id,
        assess_methods=(assess.DOPE, assess.GA341)
    )
    a.starting_model = 1
    a.ending_model = 1
    
    a.make()
    os.chdir(CURRENT_DIR)
    return "SUCCESS"

if __name__ == '__main__':
    p = Pool(cpu_count())
    print(p.map(run_modeller, allSeqs))
