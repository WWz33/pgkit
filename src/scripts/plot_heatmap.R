#!/usr/bin/env Rscript
# =============================================================================
# PAV Heatmap with ComplexHeatmap (APAVplot style)
# =============================================================================
# Features:
#   - Block splitting by category
#   - Row annotation: category + presence count
#   - Column annotation: population (optional) + presence count
#   - Raster rendering for large matrices
#
# Usage: Rscript plot_heatmap.R <pav.tsv> <out_prefix> [freq.tsv] [pop.tsv]
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

# Load packages
if (!require("ComplexHeatmap", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("ComplexHeatmap")
}
suppressPackageStartupMessages(library(ComplexHeatmap))
library(grid)

# Read PAV matrix
mat <- read.delim(pav_file, row.names = 1)
mat <- as.matrix(mat)

# Order by presence frequency
row_order <- order(rowMeans(mat), decreasing = TRUE)
mat <- mat[row_order, ]

n_genes <- nrow(mat)
n_samples <- ncol(mat)

# Determine raster for large matrices
use_raster <- n_genes * n_samples > 100000

# Colors
cat_colors <- c(
  core = "#F8766D",
  soft_core = "#7CAE00",
  dispensable = "#00BFC4",
  private = "#C77CFF"
)
pav_colors <- c("absence" = "#F0F0F0", "presence" = "#2166AC")

# ============================================================
# Row annotation
# ============================================================
df_freq <- NULL
anno_left <- NULL

if (!is.null(freq_file)) {
  df_freq <- read.delim(freq_file, row.names = 1)
  df_freq <- df_freq[rownames(mat), "Category", drop = FALSE]
  
  anno_left <- ComplexHeatmap::rowAnnotation(
    Category = ComplexHeatmap::anno_simple(
      df_freq$Category,
      col = cat_colors,
      border = FALSE,
      width = grid::unit(3, "mm")
    ),
    PN = ComplexHeatmap::anno_barplot(
      rowSums(mat),
      border = FALSE,
      width = grid::unit(10, "mm"),
      gp = grid::gpar(fill = "#2166AC", col = NA),
      axis_param = list(side = "left", gp = grid::gpar(fontsize = 7, fontface = "bold"))
    ),
    annotation_label = c("Category", "Presence\nCount"),
    annotation_name_side = "bottom",
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold")
  )
} else {
  anno_left <- ComplexHeatmap::rowAnnotation(
    PN = ComplexHeatmap::anno_barplot(
      rowSums(mat),
      border = FALSE,
      width = grid::unit(10, "mm"),
      gp = grid::gpar(fill = "#2166AC", col = NA),
      axis_param = list(side = "left", gp = grid::gpar(fontsize = 7, fontface = "bold"))
    ),
    annotation_label = "Presence\nCount",
    annotation_name_side = "bottom",
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold")
  )
}

# ============================================================
# Column annotation (top)
# ============================================================
anno_top <- NULL
pop_colors <- NULL

if (!is.null(pop_file)) {
  df_pop <- read.delim(pop_file, header = FALSE, col.names = c("species", "group"))
  rownames(df_pop) <- df_pop$species
  df_pop <- df_pop[colnames(mat), "group", drop = FALSE]
  colnames(df_pop) <- "Population"
  
  pops <- unique(df_pop$Population)
  pop_colors <- rainbow(length(pops), s = 0.6, v = 0.85)
  names(pop_colors) <- pops
  
  anno_top <- ComplexHeatmap::HeatmapAnnotation(
    Population = ComplexHeatmap::anno_simple(
      df_pop$Population,
      col = pop_colors,
      border = FALSE,
      height = grid::unit(3, "mm")
    ),
    PN = ComplexHeatmap::anno_barplot(
      colSums(mat),
      border = FALSE,
      height = grid::unit(10, "mm"),
      gp = grid::gpar(fill = "#E41A1C", col = NA),
      axis_param = list(side = "top", gp = grid::gpar(fontsize = 7, fontface = "bold"))
    ),
    annotation_label = c("Population", "Presence\nCount"),
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold")
  )
} else {
  anno_top <- ComplexHeatmap::HeatmapAnnotation(
    PN = ComplexHeatmap::anno_barplot(
      colSums(mat),
      border = FALSE,
      height = grid::unit(10, "mm"),
      gp = grid::gpar(fill = "#E41A1C", col = NA),
      axis_param = list(side = "top", gp = grid::gpar(fontsize = 7, fontface = "bold"))
    ),
    annotation_label = "Presence\nCount",
    annotation_name_side = "right",
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold")
  )
}

# ============================================================
# Split by category
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
# Legend
# ============================================================
lg <- list(
  ComplexHeatmap::Legend(
    labels = c("absence", "presence"),
    legend_gp = grid::gpar(fill = c("#F0F0F0", "#2166AC")),
    title = "PAV",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9, fontface = "bold")
  ),
  ComplexHeatmap::Legend(
    labels = names(cat_colors),
    legend_gp = grid::gpar(fill = cat_colors),
    title = "Category",
    title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9, fontface = "bold")
  )
)

if (!is.null(pop_colors)) {
  lg <- c(lg, list(
    ComplexHeatmap::Legend(
      labels = names(pop_colors),
      legend_gp = grid::gpar(fill = pop_colors),
      title = "Population",
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9, fontface = "bold")
    )
  ))
}

# ============================================================
# Draw heatmap
# ============================================================
ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "main",
  use_raster = use_raster,
  col = as.vector(pav_colors),
  show_heatmap_legend = FALSE,
  
  # Row settings
  show_row_names = FALSE,
  cluster_rows = FALSE,
  row_split = row_split,
  row_title = row_title,
  row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 0,
  row_gap = grid::unit(2, "mm"),
  row_dend_side = "left",
  
  # Column settings
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
  cluster_columns = FALSE,
  column_title_rot = 0,
  column_title_side = "top",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  
  # Annotations
  left_annotation = anno_left,
  top_annotation = anno_top,
  
  # Border
  border = TRUE,
  
  # Title
  column_title = "Gene Presence/Absence Heatmap"
)

# Save PDF
pdf(paste0(out_prefix, ".heatmap.pdf"), width = fig_width, height = fig_height)
ComplexHeatmap::draw(ht,
  heatmap_legend_list = lg,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  auto_adjust = FALSE
)
dev.off()

# Save PNG
png(paste0(out_prefix, ".heatmap.png"), width = fig_width, height = fig_height,
    units = "in", res = 300)
ComplexHeatmap::draw(ht,
  heatmap_legend_list = lg,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  auto_adjust = FALSE
)
dev.off()

cat("Saved:", paste0(out_prefix, ".heatmap.pdf/png"), "\n")
