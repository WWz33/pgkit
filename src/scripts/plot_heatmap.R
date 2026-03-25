#!/usr/bin/env Rscript
# =============================================================================
# 热图 - PAV矩阵可视化 (pheatmap)
# =============================================================================
# 用法: Rscript plot_heatmap.R <pav_matrix.tsv> <output_prefix> [frequency_table.tsv] [format]
#        format: png, pdf, svg (默认 png)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript plot_heatmap.R <pav_matrix.tsv> <output_prefix> [frequency_table.tsv] [format]\n")
  quit(status = 1)
}

pav_file <- args[1]
out_prefix <- args[2]
freq_file <- NULL
fmt <- "png"

if (length(args) >= 3 && args[3] != "NULL" && file.exists(args[3])) {
  freq_file <- args[3]
  if (length(args) >= 4) fmt <- args[4]
} else if (length(args) >= 3 && args[3] != "NULL") {
  fmt <- args[3]
}

if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}
library(pheatmap)

# 读取PAV矩阵
mat <- read.delim(pav_file, row.names = 1)
mat <- as.matrix(mat)

# 按存在频率排序(从高到低)
row_order <- order(rowMeans(mat), decreasing = TRUE)
mat <- mat[row_order, ]

n_genes <- nrow(mat)
n_samples <- ncol(mat)

# 行注释
annotation_row <- NULL
annotation_colors <- NULL

if (!is.null(freq_file)) {
  df_freq <- read.delim(freq_file, row.names = 1)
  df_freq <- df_freq[rownames(mat), "Category", drop = FALSE]
  annotation_row <- df_freq
  annotation_colors <- list(
    Category = c(
      core = "#F8766D",
      soft_core = "#7CAE00",
      dispensable = "#00BFC4",
      private = "#C77CFF"
    )
  )
}

fig_height <- max(8, n_genes * 0.015 + 4)
fig_width <- max(10, n_samples * 0.3 + 4)

# 生成热图
out_file <- paste0(out_prefix, ".heatmap.", fmt)

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
  annotation_colors = annotation_colors,
  filename = out_file,
  width = fig_width,
  height = fig_height
)

cat("Saved:", out_file, "\n")
