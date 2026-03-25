#!/usr/bin/env Rscript
# =============================================================================
# Gene Family Distribution Histogram + Ring Chart (based on APAVplot)
# =============================================================================
# Features:
#   - Histogram showing number of orthogroups by presence count
#   - Ring chart showing core/soft-core/dispensable/private proportions
#   - Combined visualization
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

# Read frequency table
df <- read.delim(freq_file)

# Category colors
cat_colors <- c(
  "core" = "#F8766D",
  "soft_core" = "#7CAE00",
  "dispensable" = "#00BFC4",
  "private" = "#C77CFF"
)

# Histogram: distribution by species count
p_hist <- ggplot(df, aes(x = Species_Count, fill = Category)) +
  geom_histogram(binwidth = 1, alpha = 0.8, color = "white") +
  scale_fill_manual(values = cat_colors) +
  labs(x = "Number of Species", y = "Number of Orthogroups",
       title = "Gene Family Distribution by Species Count") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10, color = "black", face = "bold"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

# Ring chart data
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

# Save histogram
ggsave(paste0(out_prefix, ".hist.pdf"), p_hist, width = 10, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".hist.png"), p_hist, width = 10, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".hist.pdf/png"), "\n")

# Save ring chart
ggsave(paste0(out_prefix, ".ring.pdf"), p_ring, width = 6, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".ring.png"), p_ring, width = 6, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".ring.pdf/png"), "\n")

# Combined plot (if patchwork available)
if (require("patchwork", quietly = TRUE)) {
  p_combined <- p_hist + p_ring + 
    plot_layout(widths = c(2, 1)) +
    plot_annotation(title = "Pan-genome Gene Family Analysis",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))
  ggsave(paste0(out_prefix, ".combined.pdf"), p_combined, width = 14, height = 6, dpi = 300)
  ggsave(paste0(out_prefix, ".combined.png"), p_combined, width = 14, height = 6, dpi = 300)
  cat("Saved:", paste0(out_prefix, ".combined.pdf/png"), "\n")
}
