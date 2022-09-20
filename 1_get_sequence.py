"""
Collecting sequence data and clustering from Smith et al. 2004 paper.
URL: https://science.sciencemag.org/content/305/5682/371.
"""

import os
import re
import json
from collections import defaultdict

from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from pyutils import logging, DATA_DIR, HA1_ALIGNED_FILE
from pyutils.entrez_data import PMID, ADDITIONAL, REGION_NAMES
from pyutils.miscellaneous import run_cmd


HA1_REFERENCE_FILE = os.path.join(DATA_DIR, "HA1_reference.fasta")

# Internet access to NCBI entrez service is required
Entrez.email = "youremail@domain.com"

logging.info("Use PMID to retrieve the ids of associated sequences")
# Should be only one record
(res,) = Entrez.read(Entrez.elink(linkname="pubmed_nuccore", id=PMID, idtype="acc"))
# Should be only one record
(res,) = res["LinkSetDb"]
# Add missing ids
IdList = [*[acc["Id"] for acc in res["Link"]], *ADDITIONAL]

logging.info("Download the sequences from Smith's paper")
# Download and parse in Genbank format
handle = Entrez.efetch(db="nuccore", id=IdList, rettype="gb", retmode="text")
records = [record for record in SeqIO.parse(handle, "gb")]
# Read in reference sequence
refSeq = SeqIO.read(HA1_REFERENCE_FILE, "fasta")
outSeqs = [refSeq]

logging.info("Map strain names and names in Genbank record")
seqname2ac = defaultdict(set)
for record in records:
    organism = None
    organism2 = None
    seq = None
    for feature in record.features:
        if feature.type == "source":
            # Get the organism name
            (organism,) = feature.qualifiers["organism"]
            m = re.search(r"/[A-Za-z ]+/[A-Za-z0-9]+/[0-9]+", organism)
            organism = m.group(0)
            _, region, name, year = organism.split("/")
            m = re.search(r"[0-9]+", name)
            name = m.group(0)
            # The renamed organism for mapping cluster name
            organism = (REGION_NAMES[region], name, year[-2:])
            organism2 = organism
            # Get the isolate/strain name of the virus
            if "strain" in feature.qualifiers:
                (strain,) = feature.qualifiers["strain"]
            elif "isolate" in feature.qualifiers:
                (strain,) = feature.qualifiers["isolate"]
            # Correct the year if possible
            year2 = strain.split("/")[-1]
            if year2:
                year2 = year2.split('(')[0].strip()
                organism2 = (REGION_NAMES[region], name, year2[-2:])
        if feature.type == "CDS":
            (seq,) = feature.qualifiers["translation"]
    # Map the renamed organism with accession id
    seqname2ac[organism].add(record.id)
    seqname2ac[organism2].add(record.id)
    outSeqs.append(SeqRecord(seq=Seq(seq), id=record.id, description=""))

sequences_filename = "sequences.fasta"
sequences_file = os.path.join(DATA_DIR, sequences_filename)
SeqIO.write(outSeqs, sequences_file, "fasta")

logging.info("Parse the copied cluster info")
# The file was copied and pasted
clusters = {}
with open(os.path.join(DATA_DIR, "metadata_copied.txt")) as f:
    for row in f:
        row = row.split(" ")
        seqname = row[1].split("/")
        m = re.search(r"[0-9]+", seqname[1])
        seqname[1] = m.group(0)
        if seqname[0] == "OV":
            seqname[0] = "MA"
        seqname = tuple(seqname)
        for ac in seqname2ac[seqname]:
            clusters[ac] = row[0]
with open(os.path.join(DATA_DIR, "metadata.json"), "w") as f:
    json.dump(clusters, f, indent=4)

assert all(record.id in clusters for record in records)

# Align sequences with muscle
aligned_filename = "aligned.fasta"
run_cmd("muscle", "-align", sequences_filename, "-output", aligned_filename)
aligned_file = os.path.join(DATA_DIR, "aligned.fasta")

logging.info("Remove the gaps from reference sequence for the site numbering")
# Read in the output MSA
seqs = AlignIO.read(aligned_file, "fasta")

refSeq_index = None
for n, record in enumerate(seqs):
    if record.id == refSeq.id:
        refSeq_index = n
        break

assert refSeq_index is not None

# Indices for ungapped reference sequences
prev_gap_index = -1
prev_aa_index = -1
start_index = None
end_index = None
reference_indexes = []
for index, aa in enumerate(seqs[refSeq_index]):
    if aa == "-":
        if index != prev_gap_index + 1:
            print(index, aa)
            reference_indexes.append((start_index, end_index))
        prev_gap_index = index
    else:
        if index != prev_aa_index + 1:
            start_index = index
        end_index = index + 1
        prev_aa_index = index

non_gap_seqs = None
for start_index, end_index in reference_indexes:
    if non_gap_seqs:
        non_gap_seqs += seqs[:, start_index:end_index]
    else:
        non_gap_seqs = seqs[:, start_index:end_index]

out_aligned_seqs = []
for n, record in enumerate(non_gap_seqs):
    if n != refSeq_index:
        out_aligned_seqs.append(record)

AlignIO.write(
    MultipleSeqAlignment(out_aligned_seqs),
    HA1_ALIGNED_FILE,
    "fasta"
)
