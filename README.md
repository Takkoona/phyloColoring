# Application of sitePath

The application of R package `sitePath` on SARS-CoV-2 and H3N2.

## Dependencies to run this project

Python packages: [`jupyter-notebook`](https://jupyter.readthedocs.io/en/latest/install.html), [`biopython`](https://biopython.org/wiki/Download), [`pandas`](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html). Here is one of the ways to install them.

```
pip install notebook biopython pandas matplotlib
```

R packages: [`VennDiagram`](https://cran.r-project.org/package=VennDiagram), [`gridExtra`](https://cran.r-project.org/package=gridExtra), [`svglite`](https://cran.r-project.org/package=svglite), [`jsonlite`](https://cran.r-project.org/package=jsonlite), [`IRkernel`](https://cran.r-project.org/package=IRkernel). Here is one of the ways to install them.

```r
install.packages(c("VennDiagram", "gridExtra", "svglite", "jsonlite", "IRkernel"))
```

R Bioconductor package: [`sitePath`](https://bioconductor.org/packages/sitePath/), [`ggtree`](https://bioconductor.org/packages/ggtree/), [`trackViewer`](https://bioconductor.org/packages/trackViewer/). Here is one of the ways to install them.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("sitePath", "ggtree", "trackViewer"))
```

External packages: [`muscle`]
1. Download the executable from http://www.drive5.com/muscle/ 
2. Rename the executable to "muscle" 
3. Add its location (`path/to/muscle`) to your `PATH`

## Install dependencies via bioconda on Linux & Mac

If you use [`bioconda`](https://bioconda.github.io/user/install.html#set-up-channels) on Linux or Mac, then you can install the dependencies by 

```bash
conda install notebook pandas biopython matplotlib \
    r-venndiagram r-gridextra r-svglite r-jsonlite r-irkernel \
    bioconductor-sitepath bioconductor-trackviewer \
    muscle
```
