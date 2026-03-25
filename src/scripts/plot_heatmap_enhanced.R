#!/usr/bin/env Rscript
# =============================================================================
# Enhanced PAV Heatmap (based on APAVplot concepts)
# =============================================================================
# Features:
#   - Block splitting by category (core/soft-core/dispensable/private)
#   - Row/Column bar annotations
#   - Presence number statistics
#   - Customizable colors and sizes
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_heatmap_enhanced.R <pav_matrix.tsv> <output_prefix> [frequency_table.tsv]\n")
  quit(status = 1)
}

pav_file <- args[1]
out_prefix <- args[2]
freq_file <- if (length(args) >= 3) args[3] else NULL

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

# Read frequency table for annotation
df_anno <- NULL
if (!is.null(freq_file) && file.exists(freq_file)) {
  df_freq <- read.delim(freq_file, row.names = 1)
  df_anno <- data.frame(Category = df_freq[rownames(mat), "Category"])
  rownames(df_anno) <- rownames(mat)
}

# Category colors
type_colors <- c(
  "core" = "#F8766D",
  "soft_core" = "#7CAE00",
  "dispensable" = "#00BFC4",
  "private" = "#C77CFF"
)

# PAV colors
pav_colors <- c("absence" = "#F0F0F0", "presence" = "#2166AC")

# Calculate presence number per gene
gene_presence <- rowSums(mat)

# Calculate presence number per sample
sample_presence <- colSums(mat)

n_genes <- nrow(mat)
n_samples <- ncol(mat)

# Row annotation: presence number barplot
row_anno <- rowAnnotation(
  "Presence\nCount" = anno_barplot(
    gene_presence,
    border = FALSE,
    width = unit(8, "mm"),
    gp = gpar(fill = "#2166AC", col = NA),
    axis_param = list(side = "top", gp = gpar(fontsize = 7))
  ),
  annotation_label = "Presence\nCount",
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

# Column annotation: presence number barplot
col_anno <- HeatmapAnnotation(
  "Presence\nCount" = anno_barplot(
    sample_presence,
    border = FALSE,
    height = unit(8, "mm"),
    gp = gpar(fill = "#E41A1C", col = NA),
    axis_param = list(side = "right", gp = gpar(fontsize = 7))
  ),
  annotation_label = "Presence\nCount",
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

# Category annotation (if available)
row_anno_cat <- NULL
if (!is.null(df_anno)) {
  cat_colors <- type_colors[unique(df_anno$Category)]
  row_anno_cat <- rowAnnotation(
    Category = df_anno$Category,
    col = list(Category = cat_colors),
    show_legend = TRUE,
    annotation_label = "Category",
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
  )
}

# Build heatmap
fig_height <- max(8, n_genes * 0.012 + 4)
fig_width <- max(10, n_samples * 0.25 + 5)

# Generate PDF
pdf(paste0(out_prefix, ".heatmap_enhanced.pdf"), width = fig_width, height = fig_height)

ht <- Heatmap(
  mat,
  name = "PAV",
  col = pav_colors,
  show_heatmap_legend = TRUE,
  
  # Row settings
  show_row_names = FALSE,
  cluster_rows = FALSE,
  row_dend_side = "left",
  
  # Column settings
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 6),
  cluster_columns = FALSE,
  
  # Split by category (if available)
  row_split = if (!is.null(df_anno)) df_anno$Category else NULL,
  row_title = if (!is.null(df_anno)) c("Core", "Soft-core", "Dispensable", "Private") else NULL,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 0,
  
  # Gap between blocks
  row_gap = unit(2, "mm"),
  
  # Annotations
  left_annotation = row_anno,
  right_annotation = row_anno_cat,
  top_annotation = col_anno,
  
  # Border
  border = TRUE,
  
  # Title
  column_title = "Gene Presence/Absence Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

dev.off()
cat("Saved:", paste0(out_prefix, ".heatmap_enhanced.pdf"), "\n")

# Generate PNG
png(paste0(out_prefix, ".heatmap_enhanced.png"), width = fig_width, height = fig_height, 
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Saved:", paste0(out_prefix, ".heatmap_enhanced.png"), "\n")
