#!/usr/bin/env Rscript
# =============================================================================
# Pie + Histogram Combined (APAVplot style)
# =============================================================================
# Features:
#   - Histogram showing orthogroup distribution by species count
#   - Ring chart showing category proportions
#   - Combined in one figure
#
# Usage: Rscript plot_hist_ring.R <frequency_table.tsv> <output_prefix>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_hist_ring.R <frequency_table.tsv> <output_prefix>\n")
  quit(status = 1)
}

freq_file <- args[1]
out_prefix <- args[2]

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# Read frequency table
df <- read.delim(freq_file)

# Category colors
cat_colors <- c(
  core = "#F8766D",
  soft_core = "#7CAE00",
  dispensable = "#00BFC4",
  private = "#C77CFF"
)

# ============================================================
# Histogram: distribution by species count
# ============================================================
df_hist <- df
df_hist$Category <- factor(df_hist$Category, levels = c("core", "soft_core", "dispensable", "private"))

p_hist <- ggplot(df_hist, aes(x = Species_Count, fill = Category)) +
  geom_histogram(binwidth = 1, alpha = 0.8, color = "white") +
  scale_fill_manual(values = cat_colors, name = "Category") +
  labs(x = "Number of Species", y = "Number of Orthogroups") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, color = "black", face = "bold"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

# ============================================================
# Ring chart data
# ============================================================
ring_data <- as.data.frame(table(df$Category))
colnames(ring_data) <- c("Category", "Count")
ring_data$Proportion <- ring_data$Count / sum(ring_data$Count) * 100
ring_data$Label <- paste0(ring_data$Category, "\n", ring_data$Count,
                          " (", round(ring_data$Proportion, 1), "%)")
ring_data$ymax <- cumsum(ring_data$Count)
ring_data$ymin <- c(0, head(ring_data$ymax, -1))

# Ring chart
p_ring <- ggplot(ring_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
  geom_rect(color = "white") +
  geom_text(aes(x = 3.5, y = (ymax + ymin) / 2, label = Label),
            size = 3, fontface = "bold") +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  scale_fill_manual(values = cat_colors) +
  theme_void() +
  theme(legend.position = "none")

# ============================================================
# Combined output using grid (no patchwork dependency)
# ============================================================
combined_file_png <- paste0(out_prefix, ".combined.png")
combined_file_pdf <- paste0(out_prefix, ".combined.pdf")

# Save combined PNG
png(combined_file_png, width = 14, height = 6, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(2, 1))))
print(p_hist, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_ring, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Save combined PDF
pdf(combined_file_pdf, width = 14, height = 6)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(2, 1))))
print(p_hist, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_ring, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

cat("Saved:", combined_file_pdf, "\n")
cat("Saved:", combined_file_png, "\n")
