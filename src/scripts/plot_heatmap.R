#!/usr/bin/env Rscript
# =============================================================================
# Enhanced PAV Heatmap (APAVplot style)
# =============================================================================
# Features:
#   - pheatmap with category annotation
#   - Pop (population) annotation support
#
# Usage: Rscript plot_heatmap.R <pav.tsv> <out_prefix> [freq.tsv] [pop.tsv]
#   pop.tsv: two columns (species, group), no header
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_heatmap.R <pav.tsv> <out_prefix> [freq.tsv] [pop.tsv]\n")
  quit(status = 1)
}

pav_file <- args[1]
out_prefix <- args[2]
freq_file <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL
pop_file <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL

if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}
library(pheatmap)

# Read PAV matrix
mat <- read.delim(pav_file, row.names = 1)
mat <- as.matrix(mat)

# Order by presence frequency
row_order <- order(rowMeans(mat), decreasing = TRUE)
mat <- mat[row_order, ]

n_genes <- nrow(mat)
n_samples <- ncol(mat)

# Row annotation (category)
annotation_row <- NULL
annotation_colors <- list()

if (!is.null(freq_file)) {
  df_freq <- read.delim(freq_file, row.names = 1)
  df_freq <- df_freq[rownames(mat), "Category", drop = FALSE]
  annotation_row <- df_freq
  annotation_colors[["Category"]] <- c(
    core = "#F8766D",
    soft_core = "#7CAE00",
    dispensable = "#00BFC4",
    private = "#C77CFF"
  )
}

# Column annotation (population)
annotation_col <- NULL

if (!is.null(pop_file)) {
  df_pop <- read.delim(pop_file, header = FALSE, col.names = c("species", "group"))
  rownames(df_pop) <- df_pop$species
  df_pop <- df_pop[colnames(mat), "group", drop = FALSE]
  colnames(df_pop) <- "Population"
  annotation_col <- df_pop
  
  pops <- unique(df_pop$Population)
  n_pops <- length(pops)
  pop_colors <- rainbow(n_pops, s = 0.6, v = 0.85)
  names(pop_colors) <- pops
  annotation_colors[["Population"]] <- pop_colors
}

# Figure size
fig_height <- max(8, n_genes * 0.012 + 4)
fig_width <- max(10, n_samples * 0.25 + 4)

# Generate PDF
pheatmap(mat,
  color = c("#2166AC", "#B2182B"),
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  breaks = seq(0, 1, length.out = 3),
  legend = TRUE,
  legend_breaks = c(0, 1),
  legend_labels = c("Absence", "Presence"),
  main = "Gene Presence/Absence Heatmap",
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = if (length(annotation_colors) > 0) annotation_colors else NULL,
  filename = paste0(out_prefix, ".heatmap.pdf"),
  width = fig_width,
  height = fig_height
)

# Generate PNG
pheatmap(mat,
  color = c("#2166AC", "#B2182B"),
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  breaks = seq(0, 1, length.out = 3),
  legend = TRUE,
  legend_breaks = c(0, 1),
  legend_labels = c("Absence", "Presence"),
  main = "Gene Presence/Absence Heatmap",
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = if (length(annotation_colors) > 0) annotation_colors else NULL,
  filename = paste0(out_prefix, ".heatmap.png"),
  width = fig_width,
  height = fig_height
)

cat("Saved:", paste0(out_prefix, ".heatmap.pdf/png"), "\n")
