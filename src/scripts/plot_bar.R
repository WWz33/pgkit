#!/usr/bin/env Rscript
# =============================================================================
# 柱状图 - 各样本基因家族组成
# =============================================================================
# 用法: Rscript plot_bar.R <gene_count_matrix.tsv> <output_prefix> [stacked] [format]
#        stacked: TRUE/FALSE (默认 FALSE)
#        format: png, pdf, svg (默认 png)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript plot_bar.R <gene_count_matrix.tsv> <output_prefix> [stacked] [format]\n")
  quit(status = 1)
}

count_file <- args[1]
out_prefix <- args[2]
stacked <- if (length(args) >= 3 && args[3] == "TRUE") TRUE else FALSE
fmt <- if (length(args) >= 4) args[4] else "png"

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# 读取基因计数矩阵
df <- read.delim(count_file)

# 转换为长格式
df_long <- data.frame(
  Sample = rep(df$Species, 4),
  Category = rep(c("Core", "Soft-core", "Dispensable", "Private"), each = nrow(df)),
  Count = c(df$Core, df$Soft.core, df$Dispensable, df$Private)
)

df_long$Category <- factor(df_long$Category,
  levels = c("Core", "Soft-core", "Dispensable", "Private")
)

# 按Core数量排序样本
sample_order <- df_long[df_long$Category == "Core", ]
sample_order <- sample_order[order(-sample_order$Count), "Sample"]
df_long$Sample <- factor(df_long$Sample, levels = sample_order)

# 绘图
position <- if (stacked) "stack" else "dodge"
fig_width <- max(12, nrow(df) * 0.4)

p <- ggplot(df_long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position, width = 0.7) +
  scale_fill_manual(values = c(
    "Core" = "#F8766D",
    "Soft-core" = "#7CAE00",
    "Dispensable" = "#00BFC4",
    "Private" = "#C77CFF"
  )) +
  labs(
    x = "Sample", y = "Number of Gene Families",
    title = "Gene Family Composition by Sample"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

out_file <- paste0(out_prefix, ".bar.", fmt)
ggsave(out_file, p, width = fig_width, height = 6, dpi = 300)
cat("Saved:", out_file, "\n")
