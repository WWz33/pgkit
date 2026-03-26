#!/usr/bin/env Rscript
# =============================================================================
# Stacked Bar Chart - Gene Family Composition by Sample
# =============================================================================
# Usage: Rscript plot_bar.R <gene_count_matrix.tsv> <output_prefix> [format]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_bar.R <gene_count_matrix.tsv> <output_prefix> [format]\n")
  quit(status = 1)
}

count_file <- args[1]
out_prefix <- args[2]
fmt <- if (length(args) >= 3) args[3] else "png"

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# Category colors
cat_colors <- c(
  Core = "#f78d85",
  Softcore = "#ffc725",
  Dispensable = "#48b6a6",
  Specific = "#8d9dc7"
)

# Read gene count matrix
df <- read.delim(count_file)

# Convert to long format
df_long <- data.frame(
  Sample = rep(df$Species, 4),
  Category = rep(c("Core", "Softcore", "Dispensable", "Specific"), each = nrow(df)),
  Count = c(df$Core, df$Softcore, df$Dispensable, df$Specific)
)

df_long$Category <- factor(df_long$Category,
  levels = c("Core", "Softcore", "Dispensable", "Specific")
)

# Sort samples by Core count
sample_order <- df_long[df_long$Category == "Core", ]
sample_order <- sample_order[order(-sample_order$Count), "Sample"]
df_long$Sample <- factor(df_long$Sample, levels = sample_order)

# Stacked bar chart
fig_width <- max(12, nrow(df) * 0.4)

p <- ggplot(df_long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = cat_colors) +
  labs(
    x = "Sample", y = "Number of Gene Families", fill = "Category",
    title = "Gene Family Composition by Sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

# Save
ggsave(paste0(out_prefix, ".bar.pdf"), p, width = fig_width, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".bar.png"), p, width = fig_width, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".bar.pdf/png"), "\n")
