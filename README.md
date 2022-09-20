# Test robustness of sitePath

## Download code and denpendency

This project uses [Bioconda](https://bioconda.github.io/) to manage package dependency.

```bash
conda env update
conda activate sitePath_assessment
```

## Prepare test data

To collect antigenic shift mutations of influenza A H3N2 from [Smith et al. 2004 paper](https://science.sciencemag.org/content/305/5682/371) and build tree.

```bash
python 1_get_sequence.py
python 2_build_tree.py
```

## Run test

Test sitePath using the antigenic shift mutations.

```bash
Rscript 3_sitePath_test.R
```
