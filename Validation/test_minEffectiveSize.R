library(parallel)
library(ape)
library(sitePath)

sites <- readLines("sites.txt")
sites <- as.integer(sites)

treeDir <- file.path("..", "Trees", "H1N1_HA")

tree <- read.tree(file.path(treeDir, "RAxML_bestTree.PROTGAMMAGTR"))
tree <- drop.tip(tree, "MK615591")
tree <- addMSA(tree, msaPath = file.path(treeDir, "aligned.fasta"), msaFormat = "fasta")

#---------------------------------------------------------------------------------------

testParam <- "minEffectiveSize"
dir.create(testParam, showWarnings = FALSE)

cl <- makeCluster(detectCores())
clusterExport(cl, c("tree", "testParam"))

simValues <- seq(0.01, 0.02, length.out = 3)
res <- parLapply(cl, simValues, function(similarity) {
    paths <- sitePath::lineagePath(tree, similarity)
    sizeValues <- seq(0, ape::Ntip(tree) / 4, by = 30)
    m <- lapply(sizeValues, function(ms) {
        mutations <- sitePath::fixationSites(paths, minEffectiveSize = ms)
        fileName <- file.path(testParam, paste0(similarity, "_", ms, ".rds"))
        saveRDS(mutations, fileName)
        return(mutations)
    })
    names(m) <- sizeValues
    return(m)
})
names(res) <- simValues
saveRDS(res, paste0(testParam, ".rds"))

stopCluster(cl)
