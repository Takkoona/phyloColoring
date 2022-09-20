library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(sitePath)
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(VennDiagram))


data_dir <- file.path("data", "H3N2_HA1_Smith2004")
output_dir <- "output"
plots_dir <- "plots"

phylo_tools <- list(
    "iqtree" = "gtr_ufboot.treefile",
    "raxml" = "RAxML_bestTree.PROTGAMMAGTR",
    "fasttree" = "FastTree.nwk",
    "megaML" = "MegaML_JTT.nwk"
)
# =======================================================================
# Read and process data

# Manually identified fixation sites
manual_fixation_sites <- jsonlite::read_json(file.path(
    data_dir,
    "manual_fixation_sites.json"
), simplifyVector = TRUE)

# The original paper wrote R102K for EN72-VI75,
# but I think it's a typo for R201K
cluster_mutations <- jsonlite::read_json(file.path(
    data_dir,
    "cluster_mutations.json"
), simplifyVector = TRUE)

antigenic_sites <- unique(unlist(lapply(
    cluster_mutations,
    function(mutations) {
        as.integer(substr(mutations, 2, nchar(mutations) - 1))
    }
)))

cluster_info <- jsonlite::read_json(file.path(data_dir, "metadata.json"))
cluster_info <- do.call(rbind, lapply(names(cluster_info), function(ac) {
    data.frame(
        "id" = ac,
        "cluster" = cluster_info[[ac]]
    )
}))
group_info <- split(cluster_info[["id"]], cluster_info[["cluster"]])
# Color for plot
group_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(group_info))
names(group_colors) <- names(group_info)
group_colors["0"] <- NA

alignment <- seqinr::read.alignment(
    file.path(data_dir, "HA1_aligned.fasta"),
    "fasta"
)
# =======================================================================

# Build target dataset
sequences <- alignment$seq
loci <- which(vapply(
    X = seq_len(nchar(sequences[[1]])),
    FUN = function(s) {
        length(unique(substr(sequences, s, s))) > 1
    },
    FUN.VALUE = logical(1)
))
conserved_sites <- setdiff(loci, manual_fixation_sites)

site_category <- data.frame(
    "site" = c(manual_fixation_sites, conserved_sites),
    "category" = c(
        rep("fixed", length(manual_fixation_sites)),
        rep("conserved", length(conserved_sites))
    ),
    "fixationSite" = c(
        rep(TRUE, length(manual_fixation_sites)),
        rep(FALSE, length(conserved_sites))
    )
)
# Sort by site position
site_category <- site_category[order(site_category[["site"]]), ]
row.names(site_category) <- NULL

n_positive <- sum(site_category[["fixationSite"]])
n_negative <- sum(!site_category[["fixationSite"]])

# =======================================================================

# iqtree test
toolname <- "iqtree"

tree <- read.tree(file.path(data_dir, phylo_tools[[toolname]]))
tree <- groupOTU(tree, group_info, group_name = "cluster")

grouplevel <- levels(attr(tree, "cluster"))
grouplevel <- grouplevel[order(sapply(
    gsub("[^0-9]", "", grouplevel),
    function(year) {
        if (as.integer(year) > 50) {
            res <- paste0(19, year)
        } else {
            res <- paste0(20, year)
        }
        as.integer(res)
    }
), decreasing = TRUE)]

attr(tree, "cluster") <- factor(attr(tree, "cluster"), grouplevel)
p <- ggtree(tree, aes(color = cluster)) +
    scale_color_manual(
        values = as.list(group_colors),
        limits = setdiff(grouplevel, "0"),
        na.translate = TRUE,
        na.value = "black",
        name = "Antigenic\ncluster"
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 3)
        )
    ) +
    theme(
        plot.margin = unit(c(5.5, 5.5, 5.5, 128), "pt"),
        legend.position = c(-0.3, 0.5)
    ) +
    theme(legend.position = "left")
p

ggsave(
    filename = file.path(
        plots_dir,
        paste0("Smith2004_tree_", toolname, ".svg")
    ),
    plot = p,
    width = 6.7,
    height = 4
)

# Test sitePath with iqtree
paths <- addMSA(tree, alignment = alignment)
mutations <- fixationSites(paths)

# Venn plot between fixation mutations and antigenic sites
venn_plot <- draw.pairwise.venn(
    area1 = length(mutations),
    area2 = length(antigenic_sites),
    cross.area = length(intersect(
        antigenic_sites,
        as.integer(allSitesName(mutations))
    )),
    category = c("Fixation", ""),
    ext.text = FALSE,
    offset = 1,
    fill = c("blue", "red"),
    alpha = c(0.1, 0.1),
    lwd = c(3, 3),
    col = c("blue", "red"),
    label.col = c("blue", "black", "red"),
    cat.col = c("blue", "red"),
    cat.just = list(c(1, -1), c(1, -3)),
    rotation.degree = 180,
    margin = 0.2
)
grid.draw(venn_plot)

svg(
    filename = file.path(plots_dir, "Smith2004_venn.svg"),
    width = 3,
    height = 3,
    bg = "transparent"
)
grid.draw(venn_plot)
invisible(dev.off())

# sitePath performance against manually identified fixation sites
pred_result <- data.frame(
    "Nmin" = integer(),
    "rate" = double(),
    "category" = character()
)

for (Nmin in seq_len(20)[-1]) {
    paths <- lineagePath(paths, Nmin)
    mutations <- fixationSites(paths)

    assess_table <- site_category
    sites <- as.integer(allSitesName(mutations))
    assess_table[["predFixed"]] <- assess_table[["site"]] %in% sites
    x <- assess_table[["fixationSite"]] + assess_table[["predFixed"]]
    senstivity <- length(which(x == 2)) / n_positive
    specificity <- length(which(x == 0)) / n_negative
    pred_result <- rbind(
        pred_result,
        data.frame(
            "Nmin" = c(Nmin, Nmin),
            "rate" = c(senstivity, specificity),
            "category" = c("Senstivity", "Specificity")
        )
    )
}
pred_result <- pred_result[order(pred_result[["Nmin"]]), ]
row.names(pred_result) <- NULL

p <- ggplot(pred_result, aes(Nmin, rate)) +
    geom_point(aes(color = category)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
    labs(x = "Nmin", y = "", color = "Metric") +
    geom_vline(
        xintercept = c(5, 10, 12, 15),
        alpha = 0.3,
        linetype = "dotted"
    ) +
    theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(color = "black", fill = "transparent"),
        panel.grid = element_blank(), # All grid lines
        legend.title = element_blank()
    )
p

ggsave(
    filename = file.path(plots_dir, "Nmin.svg"),
    plot = p,
    device = "svg",
    width = 3.25, height = 2
)

# =======================================================================

# Other tools
pred_result_2 <- data.frame(
    "software" = integer(),
    "rate" = double(),
    "category" = character()
)

for (toolname in names(phylo_tools)) {
    tree <- read.tree(file.path(data_dir, phylo_tools[[toolname]]))
    tree <- ape::root(tree, "AF201874.1")
    p <- ggtree(
        groupOTU(tree, group_info, group_name = "Cluster"),
        aes(color = Cluster)
    ) +
        scale_color_manual(
            values = as.list(group_colors),
            limits = setdiff(names(group_colors), "0"),
            na.translate = TRUE,
            na.value = "black"
        )
    print(p)
    paths <- addMSA(tree, alignment = alignment)
    mutations <- fixationSites(paths)

    assess_table <- site_category
    sites <- as.integer(allSitesName(mutations))
    assess_table[["predFixed"]] <- assess_table[["site"]] %in% sites
    x <- assess_table[["fixationSite"]] + assess_table[["predFixed"]]
    senstivity <- length(which(x == 2)) / n_positive
    specificity <- length(which(x == 0)) / n_negative
    pred_result_2 <- rbind(
        pred_result_2,
        data.frame(
            "software" = c(toolname, toolname),
            "rate" = c(senstivity, specificity),
            "category" = c("Senstivity", "Specificity")
        )
    )
}

pred_result_2$software <- factor(
    pred_result_2$software,
    levels = names(phylo_tools)
)

p <- ggplot(pred_result_2, aes(x = software, y = rate, fill = category)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    #     ylim(0, 1) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
    labs(x = "Software", y = "", fill = "Metric") +
    theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(color = "black", fill = "transparent"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

p

ggsave(
    filename = file.path(plots_dir, "software.svg"),
    plot = p,
    device = "svg",
    width = 3.25, height = 2.4
)
