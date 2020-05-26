#!/usr/bin/python
#$ -V
#$ -cwd
#$ -S /usr/bin/python
#$ -q NGS
#$ -pe mpi 1
#$ -N H1N1_HA_modeller

import os

import modeller
from modeller.automodel import automodel, assess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

PDB_DIR = "."
MODELLING_DIR = os.path.join(PDB_DIR, "Modelling")

TEMPLATE_ID = "4lxv"

env = modeller.environ()

cwd = os.getcwd()

for record in SeqIO.parse(os.path.join(PDB_DIR, "H1N1_HA.fasta"), "fasta"):
    record.description = ':'.join(["sequence", record.id, "", "", "", "", "", "", "0.00", "0.00"])
    record.name = ""
    record.seq.alphabet = generic_protein
    
    outDir = os.path.join(MODELLING_DIR, record.id)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    os.chdir(outDir)
    
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
    os.chdir(cwd)
