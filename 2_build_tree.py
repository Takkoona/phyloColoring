from pyutils import HA1_ALIGNED_FILENAME, OUTGROUP
from pyutils.miscellaneous import run_cmd


run_cmd(
    "fasttree",
    "-gamma",
    "-boot", "1000",
    "-out", "FastTree.nwk",
    HA1_ALIGNED_FILENAME
)

run_cmd(
    "iqtree",
    "-pre", "gtr_ufboot",
    "-m", "GTR20+I+G4",
    "-s", HA1_ALIGNED_FILENAME,
    "-o", OUTGROUP,
    "-bb", "1000",
    "-bnni",
    "-nt", "AUTO"
)

model_name = "PROTGAMMAGTR"
run_cmd(
    "raxmlHPC",
    "-f", "a",
    "-p", "12345",
    "-x", "12345",
    "-#", "1000",
    "-s", HA1_ALIGNED_FILENAME,
    "-n", model_name,
    "-m", model_name,
    "-o", OUTGROUP
)
