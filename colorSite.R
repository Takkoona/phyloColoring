library(seqinr)
library(ggplot2)
library(ggtree)

tree <- read.tree("./Data/H1N1/HA1/Before2009/RAxML_bestTree.H1N1") # read tree in newick format
align <- read.alignment("./Data/H1N1/HA1/Before2009/aligned.fasta", "fasta") # read alignment in fasta format

saveFile <- "./Figures/test.png" # figure saving file path
singleSite <- 156 # the site to be colored on tree

plotSingleSite <- function(tree, align, singleSite) {
  # hex code for each amino acid
  AA_COLOR <- c(
    His = "#8282D2", Arg = "#9370DB", Lys = "#145AFF",
    Ile = "#55AE3A", Phe = "#3232AA", Leu = "#0F820F",
    Trp = "#B45AB4", Ala = "#C8C8C8", Met = "#FFD700",
    Pro = "#DC9682", Val = "#2F4F2F", Asn = "#00DCDC",
    Cys = "#E6E600", Gly = "#666666", Ser = "#FF6347",
    Tyr = "#ADD8E6", Gln = "#0099CC", Thr = "#FA9600",
    Glu = "#8C1717", Asp = "#E60A0A", gap = "#000000",
    unknown = "#000000", Ile_or_Leu = "#000000",
    Asp_or_Asn = "#000000", Glu_or_Gln = "#000000"
  )
  # transform one-letter amino acid to three-letter abbreviation
  AA_FULL_NAME = c(
    h = "His", r = "Arg", k = "Lys",
    i = "Ile", f = "Phe", l = "Leu",
    w = "Trp", a = "Ala", m = "Met",
    p = "Pro", v = "Val", n = "Asn",
    c = "Cys", g = "Gly", s = "Ser",
    y = "Tyr", q = "Gln", t = "Thr",
    e = "Glu", d = "Asp", `-` = "gap",
    x = "unknown", j = 'Ile_or_Leu',
    b = 'Asp_or_Asn', z = 'Glu_or_Gln'
  )
  # split sequence into vector of letters as R can't literate string
  align <- strsplit(tolower(align$seq[pmatch(tree$tip.label, align$nam)]), "")
  # get amnio acid variation for 'singleSite' in the alignment
  siteComp <- sapply(align, "[[", singleSite)
  # map the amnio acid variation with the tree tip name
  names(siteComp) <- tree$tip.label
  # initialize an emtpy list to store grouping
  group <- list()
  # iterate and group
  for (s in names(siteComp)) {
    # group the tree tips if they're the same amino acid at 'singleSite'
    group[[siteComp[[s]]]] <- c(group[[siteComp[[s]]]], s)
  }
  # convert one-letter amino acid to three-letter abbreviation
  names(group) <- AA_FULL_NAME[names(group)]
  # get the hex color code for each amino acid
  groupCols <- AA_COLOR[names(group)]
  # use ggptree to plot and color the tree
  p <- ggtree(groupOTU(tree, group), aes(color = group)) +
    scale_color_manual(values = c('white', groupCols)) +
    theme(legend.position = "left") +
    ggtitle(singleSite) +
    guides(
      color = guide_legend(
        override.aes = list(size = 3),
        title = 'Amino acid'
      )
    )
  return(p)
}

p <- plotSingleSite(tree, align, singleSite) # specify 'tree', 'align', 'plotSingleSite' to plot tree

ggsave(saveFile, p, device = "png", width = 7, height = 10) # save the plot to 'saveFile'
