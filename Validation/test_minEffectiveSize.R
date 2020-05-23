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

simValues <- seq(0.001, 0.01, length.out = 3)

res <- list()
for (similarity in simValues) {
    paths <- sitePath::lineagePath(tree, similarity)

    cl <- makeCluster(detectCores(), outfile = "log.txt")
    clusterExport(cl, c("paths", "testParam", "similarity"))
    
    sizeValues <- seq(0, ape::Ntip(tree) / 4, by = 30)
    m <- parLapply(cl, sizeValues, function(ms) {
        fileName <- file.path(testParam, paste0(similarity, "_", ms, ".rds"))
        cat(fileName, "\n")
        mutations <- sitePath::fixationSites(paths, minEffectiveSize = ms, method = "insert")
        saveRDS(mutations, fileName)
        return(mutations)
    })
    names(m) <- sizeValues
    res[[as.character(similarity)]] <- m
    
    stopCluster(cl)
}

saveRDS(res, paste0(testParam, ".rds"))
