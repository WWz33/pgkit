#!/usr/bin/env Rscript
# =============================================================================
# PAV Heatmap with ComplexHeatmap (APAVplot style)
# =============================================================================
# Features:
#   - Block splitting by category (Core/Soft-core/Dispensable/Private)
#   - Row annotation: category + presence count bar
#   - Column annotation: population + presence count bar
#   - Population annotation from pop file
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

# Install ComplexHeatmap if not available
if (!require("ComplexHeatmap", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)
library(grid)

# Read PAV matrix
mat <- read.delim(pav_file, row.names = 1)
mat <- as.matrix(mat)

# Order by presence frequency (high to low)
row_order <- order(rowMeans(mat), decreasing = TRUE)
mat <- mat[row_order, ]

n_genes <- nrow(mat)
n_samples <- ncol(mat)

# Category data
df_freq <- NULL
if (!is.null(freq_file)) {
  df_freq <- read.delim(freq_file, row.names = 1)
  df_freq <- df_freq[rownames(mat), "Category", drop = FALSE]
}

# Category colors (APAVplot style)
cat_colors <- c(
  core = "#F8766D",
  soft_core = "#7CAE00",
  dispensable = "#00BFC4",
  private = "#C77CFF"
)

# PAV colors
pav_colors <- c("0" = "#F0F0F0", "1" = "#2166AC")

# ============================================================
# Row annotation: category + presence count
# ============================================================
row_anno_list <- list()

# Category bar (left side)
if (!is.null(df_freq)) {
  row_anno_list[["Category"]] <- anno_simple(
    df_freq$Category,
    col = cat_colors,
    border = FALSE,
    width = unit(3, "mm")
  )
}

# Presence count bar (left side)
row_anno_list[["Presence\nCount"]] <- anno_barplot(
  rowSums(mat),
  border = FALSE,
  width = unit(10, "mm"),
  gp = gpar(fill = "#2166AC", col = NA),
  axis_param = list(side = "top", gp = gpar(fontsize = 7))
)

row_anno <- do.call(rowAnnotation, row_anno_list)

# ============================================================
# Column annotation: population + sample count
# ============================================================
col_anno_list <- list()

# Population annotation (top side)
pop_colors <- NULL
if (!is.null(pop_file)) {
  df_pop <- read.delim(pop_file, header = FALSE, col.names = c("species", "group"))
  rownames(df_pop) <- df_pop$species
  df_pop <- df_pop[colnames(mat), "group", drop = FALSE]
  colnames(df_pop) <- "Population"

  pops <- unique(df_pop$Population)
  pop_colors <- rainbow(length(pops), s = 0.6, v = 0.85)
  names(pop_colors) <- pops

  col_anno_list[["Population"]] <- anno_simple(
    df_pop$Population,
    col = pop_colors,
    border = FALSE,
    height = unit(3, "mm")
  )
}

# Sample count bar (top side)
col_anno_list[["Presence\nCount"]] <- anno_barplot(
  colSums(mat),
  border = FALSE,
  height = unit(10, "mm"),
  gp = gpar(fill = "#E41A1C", col = NA),
  axis_param = list(side = "right", gp = gpar(fontsize = 7))
)

col_anno <- do.call(HeatmapAnnotation, col_anno_list)

# ============================================================
# Split by category (block splitting)
# ============================================================
row_split <- NULL
row_title <- NULL
if (!is.null(df_freq)) {
  row_split <- factor(df_freq$Category, levels = c("core", "soft_core", "dispensable", "private"))
  row_title <- c("Core", "Soft-core", "Dispensable", "Private")
}

# ============================================================
# Figure size
# ============================================================
fig_height <- max(8, n_genes * 0.012 + 4)
fig_width <- max(10, n_samples * 0.25 + 5)

# ============================================================
# Draw heatmap
# ============================================================
ht <- Heatmap(
  mat,
  name = "PAV",
  col = pav_colors,
  show_heatmap_legend = TRUE,

  # Row settings
  show_row_names = FALSE,
  cluster_rows = FALSE,
  row_split = row_split,
  row_title = row_title,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),

  # Column settings
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 6),
  cluster_columns = FALSE,

  # Annotations
  left_annotation = row_anno,
  top_annotation = col_anno,

  # Title
  column_title = "Gene Presence/Absence Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

# Save PDF
pdf(paste0(out_prefix, ".heatmap.pdf"), width = fig_width, height = fig_height)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Save PNG
png(paste0(out_prefix, ".heatmap.png"), width = fig_width, height = fig_height,
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("Saved:", paste0(out_prefix, ".heatmap.pdf/png"), "\n")
