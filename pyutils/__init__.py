import os
import logging

DATA_DIR = os.path.join("data", "H3N2_HA1_Smith2004")

HA1_ALIGNED_FILENAME = "HA1_aligned.fasta"
HA1_ALIGNED_FILE = os.path.join(DATA_DIR, HA1_ALIGNED_FILENAME)
OUTGROUP = "AF201874.1"

# Logging config
logging.basicConfig(
    format="[%(asctime)s %(process)s]: %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO
)
