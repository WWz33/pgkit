#!/usr/bin/env Rscript
# =============================================================================
# Histogram + Ring Chart Combined (APAVplot style)
# =============================================================================
# Usage: Rscript plot_hist_ring.R <frequency_table.tsv> <output_prefix>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_hist_ring.R <frequency_table.tsv> <output_prefix>\n")
  quit(status = 1)
}

freq_file <- args[1]
out_prefix <- args[2]

# Load packages
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Source palette
script_dir <- if (is.null(sys.frame(1)$ofile)) "." else dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "palette.R"))

# Category colors and levels
cat_colors <- CATEGORY_COLORS
cat_levels <- names(CATEGORY_COLORS)

# Read frequency table
df <- read.delim(freq_file)

# Category order
cat_levels <- c("Core", "Softcore", "Dispensable", "Private")
df$Category <- factor(df$Category, levels = cat_levels)

# ============================================================
# Histogram: distribution by species count
# ============================================================
p_hist <- ggplot() +
  geom_bar(data = df, aes(x = Species_Count, y = ..count.., 
                          fill = factor(Category, levels = cat_levels)),
           stat = "count", width = 0.8) +
  theme_bw() +
  scale_fill_manual(values = cat_colors) +
  labs(x = "Number of Species", y = "Count", fill = "Category") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black", face = "bold"),
    axis.title = element_text(size = 12, color = "black", face = "bold")
  )

# ============================================================
# Ring chart data
# ============================================================
ring_data <- as.data.frame(table(df$Category))
colnames(ring_data) <- c("Category", "Count")
ring_data$Category <- factor(ring_data$Category, levels = cat_levels)
ring_data <- ring_data[!is.na(ring_data$Count), ]

ring_data$x <- ring_data$Count / 2 + c(0, head(cumsum(ring_data$Count), -1))
ring_data$per <- ring_data$Count / sum(ring_data$Count) * 100
ring_data$Label <- paste0(ring_data$Category, "\n(", ring_data$Count, ", ", 
                          round(ring_data$per, 1), "%)")

# Ring chart
p_ring <- ggplot(ring_data, aes(x = x, y = 1, width = Count, height = 2)) +
  geom_tile(aes(fill = Category)) +
  geom_text(aes(y = 2.5, label = Label), size = 3, fontface = "bold") +
  coord_polar("x") +
  scale_y_continuous(limits = c(-2, 4)) +
  scale_fill_manual(values = cat_colors) +
  theme_void() +
  theme(legend.position = "none")

# ============================================================
# Combined output using grid viewport
# ============================================================
combined_png <- paste0(out_prefix, ".combined.png")
combined_pdf <- paste0(out_prefix, ".combined.pdf")

# Save PNG
png(combined_png, width = 10, height = 6, units = "in", res = 300)
grid::grid.newpage()
print(p_hist)
vp <- grid::viewport(x = 0.65, y = 0.65, width = 0.35, height = 0.5)
print(p_ring, vp = vp)
dev.off()

# Save PDF
pdf(combined_pdf, width = 10, height = 6)
grid::grid.newpage()
print(p_hist)
vp <- grid::viewport(x = 0.65, y = 0.65, width = 0.35, height = 0.5)
print(p_ring, vp = vp)
dev.off()

cat("Saved:", combined_png, "\n")
cat("Saved:", combined_pdf, "\n")
