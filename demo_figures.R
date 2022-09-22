library(sitePath)
library(ggplot2)
suppressPackageStartupMessages(library(ggtree))

output_plot_dir <- file.path("plots", "demo_figures")
tree_height_ratio <- 1.5

if (dir.exists(output_plot_dir)) {
    unlink(output_plot_dir, recursive = TRUE)
}
dir.create(output_plot_dir, showWarnings = FALSE)

plot_path_comparison <- function(paths,
                                 select1,
                                 select2,
                                 path_size = NULL,
                                 branch_size = 1) {
    if (is.null(path_size)) {
        path_size <- 6
    }
    path_nodes <- unique(unlist(paths[c(select1, select2)]))

    group <- rep(1, times = n_nodes)
    group[path_nodes] <- 0
    group <- factor(group)

    group_colors <- rep(2, times = n_nodes)
    group_colors[path_nodes] <- 1
    group_colors[unlist(paths[select1])] <- 0
    group_colors <- factor(group_colors)

    size <- rep(1, times = n_nodes)
    size[path_nodes] <- 2

    p <- ggtree(tree, aes(
        color = group_colors,
        linetype = group,
        size = size
    )) +
        scale_size(range = c(branch_size, path_size)) +
        scale_color_manual(values = c("red", "blue", "gainsboro")) +
        theme(legend.position = "none")
    p
}

plot_sub_paths <- function(paths,
                           select = NULL,
                           path_color = "blue",
                           path_size = 2,
                           branch_size = 1) {
    tree <- ape::as.phylo(paths)
    if (is.null(select)) {
        path_nodes <- unique(unlist(paths))
    } else {
        path_nodes <- unique(unlist(paths[select]))
    }
    group <- rep(1, times = n_nodes)
    group[path_nodes] <- 0
    group <- factor(group)
    size <- rep(1, times = n_nodes)
    size[path_nodes] <- 2
    p <- ggtree(tree, aes(color = group, linetype = group, size = size)) +
        scale_size(range = c(branch_size, path_size)) +
        scale_color_manual(values = c(path_color, "gainsboro")) +
        theme(legend.position = "none")
    p
}

# =======================================================================

data(zikv_align)
data(zikv_tree)

zikv_paths <- addMSA(zikv_tree, alignment = zikv_align)
tree <- as.phylo(zikv_paths)
n_nodes <- ape::Nnode(tree, internal.only = FALSE)

pdf(file = file.path(output_plot_dir, "sneakPeek.pdf"))
sp <- sneakPeek(zikv_paths)
invisible(dev.off())

zikv_paths <- lineagePath(zikv_paths, 16)
zikv_min_entropy <- sitesMinEntropy(zikv_paths)
zikv_para_sites <- parallelSites(zikv_min_entropy, minSNP = 1)
zikv_fixed_sites <- fixationSites(zikv_min_entropy)

intersect(allSitesName(zikv_fixed_sites), allSitesName(zikv_para_sites))

# =======================================================================

p_entropy <- plotSingleSite(zikv_min_entropy, site = 139)
p_paths_mutations <- plotSingleSite(zikv_fixed_sites, site = 139)

ggsave(
    filename = file.path(output_plot_dir, "ZIKV_example.svg"),
    plot = gridExtra::arrangeGrob(p_entropy, p_paths_mutations, ncol = 2),
    device = "svg",
    width = 10,
    height = 6
)

p_mutations <- plotSingleSite(zikv_para_sites, site = 1118) +
    ggtitle(label = NULL, subtitle = NULL) +
    theme(legend.position = "none")

ggsave(
    filename = file.path(output_plot_dir, "mutations.svg"),
    plot = p_mutations,
    device = "svg",
    width = 3.5,
    height = 3.5 * tree_height_ratio
)

# =======================================================================

select1 <- 2

p <- plot_path_comparison(
    zikv_paths,
    select1,
    seq_along(zikv_paths)[-select1],
    path_size = 1.5,
    branch_size = 0.5
)
p

ggsave(
    filename = file.path(output_plot_dir, "pathComparison0.svg"),
    plot = p,
    device = "svg",
    width = 3.5,
    height = 3.5 * tree_height_ratio
)

n <- 0

for (select2 in seq_along(zikv_paths)[-select1]) {
    n <- n + 1
    p <- plot_path_comparison(zikv_paths, select1, select2, 1, 0.05)
    ggsave(
        filename = file.path(
            output_plot_dir,
            paste0("pathComparison", n, ".svg")
        ),
        plot = p,
        device = "svg",
        width = 1,
        height = 1 * tree_height_ratio
    )
}
