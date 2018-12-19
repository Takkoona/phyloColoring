library(seqinr)
library(ggplot2)
library(ggtree)

groupDir <- "./Data/H1N1/HA1/Before2009/group/" # the directory having the group information.
tree <- read.tree("./Data/H1N1/HA1/Before2009/RAxML_bestTree.H1N1") # read tree in newick format
saveFile <- "./Figures/test_grouping.png" # figure saving file path

plotGrouping <- function(tree, groupDir) {
  # initialize an empty list to store grouping
  group <- list()
  # iterate and group
  for (g in list.files(groupDir)) {
    # group was given in fasta files (same group, same file)
    # file names are the group names
    seqs <- suppressWarnings(read.alignment(file.path(groupDir, g), "fasta"))
    # extract sequence namnes and put them under the group name
    group[[g]] <- seqs$nam
  }
  # use ggtree to plot and color the tree
  p <- ggtree(groupOTU(tree, group), aes(color = group)) +
    geom_tiplab(size = 0.25) +
    theme(legend.position = "left") +
    ggtitle("Antigenic grouping") +
    guides(
      color = guide_legend(
        override.aes = list(size = 3),
        title = 'Group'
      )
    )
  return(p)
}

p <- plotGrouping(tree, groupDir) # specift 'tree' and 'groupDir' to plot tree

ggsave(saveFile, p, device = "png", width = 7, height = 10) # save the plot to 'saveFile'
