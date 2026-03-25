#!/usr/bin/env Rscript
# =============================================================================
# 饼图 - 各类别基因家族比例
# =============================================================================
# 用法: Rscript plot_pie.R <frequency_table.tsv> <output_prefix> [format]
#        format: png, pdf, svg (默认 png)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript plot_pie.R <frequency_table.tsv> <output_prefix> [format]\n")
  quit(status = 1)
}

freq_file <- args[1]
out_prefix <- args[2]
fmt <- if (length(args) >= 3) args[3] else "png"

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# 读取频率表
df <- read.delim(freq_file)

# 统计各类别
cat_counts <- table(df$Category)
cat_df <- data.frame(
  Category = names(cat_counts),
  Count = as.integer(cat_counts)
)
cat_df$Proportion <- cat_df$Count / sum(cat_df$Count)
cat_df$Label <- paste0(cat_df$Category, "\n", cat_df$Count, " (", round(cat_df$Proportion * 100, 1), "%)")

# 绘图
p <- ggplot(cat_df, aes(x = "", y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_manual(values = c(
    "core" = "#F8766D",
    "soft_core" = "#7CAE00",
    "dispensable" = "#00BFC4",
    "private" = "#C77CFF"
  )) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
  ggtitle("Gene Family Category Distribution") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

out_file <- paste0(out_prefix, ".pie.", fmt)
ggsave(out_file, p, width = 8, height = 6, dpi = 300)
cat("Saved:", out_file, "\n")
