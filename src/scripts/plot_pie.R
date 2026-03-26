#!/usr/bin/env Rscript
# =============================================================================
# Pie Chart - Gene Family Category Proportions
# =============================================================================
# Usage: Rscript plot_pie.R <frequency_table.tsv> <output_prefix> [format]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_pie.R <frequency_table.tsv> <output_prefix> [format]\n")
  quit(status = 1)
}

freq_file <- args[1]
out_prefix <- args[2]
fmt <- if (length(args) >= 3) args[3] else "png"

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Default category colors (can be overridden by palette.R)
CATEGORY_COLORS <- c(
  Core = "#f78d85",
  Softcore = "#ffc725",
  Dispensable = "#48b6a6",
  Private = "#8d9dc7"
)

# Source palette if available
tryCatch(source(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "palette.R")), error = function(e) NULL)

# Read frequency table
df <- read.delim(freq_file)

# Count categories
cat_counts <- table(df$Category)
cat_df <- data.frame(
  Category = names(cat_counts),
  Count = as.integer(cat_counts)
)
cat_df$Proportion <- cat_df$Count / sum(cat_df$Count)
cat_df$Label <- paste0(cat_df$Category, "\n", cat_df$Count, " (", round(cat_df$Proportion * 100, 1), "%)")

# Get colors for categories present in data
cat_colors <- CATEGORY_COLORS[names(cat_counts)]

# Plot
p <- ggplot(cat_df, aes(x = "", y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_manual(values = cat_colors) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
  ggtitle("Gene Family Category Distribution") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(paste0(out_prefix, ".pie.", fmt), p, width = 8, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".pie.", fmt), "\n")
