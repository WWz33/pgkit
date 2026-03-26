#!/usr/bin/env Rscript
# =============================================================================
# Half-Violin + Jitter Plot (based on APAVplot)
# =============================================================================
# Features:
#   - Half violin showing density distribution
#   - Jitter points for individual samples
#   - Clean and compact visualization
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_halfviolin.R <gene_count_matrix.tsv> <output_prefix>\n")
  quit(status = 1)
}

count_file <- args[1]
out_prefix <- args[2]

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Read data
df <- read.delim(count_file)

# Calculate total gene count per sample
sample_data <- data.frame(
  sample = df$Species,
  total = df$Total,
  stringsAsFactors = FALSE
)

# Density data for half violin
den_data <- data.frame(
  loc = density(sample_data$total)$x,
  den = density(sample_data$total)$y
)
den_data$den <- den_data$den / max(den_data$den) / 2
den_data <- subset(den_data, loc >= min(sample_data$total) & loc <= max(sample_data$total))
den_data <- rbind(den_data, c(max(sample_data$total), 0), c(min(sample_data$total), 0))

# Half violin + jitter
p <- ggplot() +
  geom_polygon(data = den_data, aes(x = -den, y = loc), fill = "#7e9bc0", alpha = 0.7) +
  geom_jitter(data = sample_data, aes(x = 0.25, y = total), 
              width = 0.125, size = 2, alpha = 0.6, color = "#E41A1C") +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  labs(y = "Total Gene Count", title = "Gene Count Distribution") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

ggsave(paste0(out_prefix, ".halfviolin.pdf"), p, width = 4, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".halfviolin.png"), p, width = 4, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".halfviolin.pdf/png"), "\n")
