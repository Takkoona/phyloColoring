#!/bin/Rscript

library(ape)
library(sitePath)

setwd("./Data/20200901/")

#Set multi-core.
options(cl.cores = 8) #The number equal the core for calculation.

#Load tree file. 
tree <- read.tree('Tree.nwk')

#Load metadata for treeID conversion and convert them. Optional, when the ID in tree and MSA file are not the same.
meta_tree <- as.data.frame(read.csv('meta_file.csv', header = T))
label_order <- as.data.frame(tree$tip.label)
colnames(label_order) <- 'Strain' #Set the column name as this string. 
new_tab <- merge(label_order,meta_tree,all.x =TRUE,sort=FALSE)
tree$tip.label <- new_tab$gisaid_epi_isl

#Reduced the tips that have no matched sequence. Optional, when the seqs in MSA file are less than the tips in tree.
reduced_list <- read.csv('Trim_ID.csv', header = F)
for (n in reduced_list){
    cat(n, "\n")
    tree <- drop.tip(tree,n)
}

cat("Calculating lineage paths...\n")
#Assign the tree with MSA.
paths <- addMSA(tree, 'aligned.fasta', "fasta")

cat("Entropy minimization...\n")
#Calculate minEntropy.
minEntropy <- sitesMinEntropy(paths)

cat("Save as RDS file\n")
#Save the calculated minEntrophy.
saveRDS(minEntropy, file.path("..", '2020-09-01.rds'))
