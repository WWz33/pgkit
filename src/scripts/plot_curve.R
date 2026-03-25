#!/usr/bin/env Rscript
# =============================================================================
# 饱和曲线 - 核心/泛基因家族
# =============================================================================
# 用法: Rscript plot_curve.R <saturation_curve.tsv> <output_prefix>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("用法: Rscript plot_curve.R <saturation_curve.tsv> <output_prefix>\n")
  quit(status = 1)
}

data_file <- args[1]
out_prefix <- args[2]

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# 读取饱和曲线数据
df <- read.delim(data_file)

# 转换为长格式
df_long <- data.frame(
  n_accessions = rep(df$n_accessions, 2),
  type = rep(c("Core", "Pan"), each = nrow(df)),
  count = c(df$core, df$pan)
)

# 绘图
p <- ggplot(df_long, aes(x = n_accessions, y = count, color = type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Core" = "#2E86AB", "Pan" = "#A23B72")) +
  labs(
    x = "Number of Accessions",
    y = "Number of Gene Families",
    color = "Type",
    title = "Core/Pan Gene Family Saturation Curve"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(paste0(out_prefix, ".curve.pdf"), p, width = 10, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".curve.png"), p, width = 10, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".curve.pdf/png"), "\n")
