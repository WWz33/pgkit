#!/usr/bin/env Rscript
# =============================================================================
# Enhanced Stacked Barplot with Dendrogram (based on APAVplot concepts)
# =============================================================================
# Features:
#   - Stacked bar by category
#   - Hierarchical clustering dendrogram
#   - Relative or absolute values
#   - Customizable colors and sizes
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_stackbar_enhanced.R <gene_count_matrix.tsv> <output_prefix> [relative]\n")
  quit(status = 1)
}

count_file <- args[1]
out_prefix <- args[2]
show_relative <- if (length(args) >= 3 && args[3] == "TRUE") TRUE else FALSE

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!require("ggdendro", quietly = TRUE)) {
  install.packages("ggdendro", repos = "https://cloud.r-project.org")
}
library(ggplot2)
library(ggdendro)

# Read data
df <- read.delim(count_file)

# Categories
categories <- c("Core", "Soft-core", "Dispensable", "Private")
cat_colors <- c("Core" = "#F8766D", "Soft-core" = "#7CAE00", 
                "Dispensable" = "#00BFC4", "Private" = "#C77CFF")

# Prepare data
p_data <- df[, c("Core", "Soft.core", "Dispensable", "Private")]
colnames(p_data) <- categories
rownames(p_data) <- df$Species

if (show_relative) {
  p_data <- as.data.frame(t(apply(p_data, 1, function(x) x / sum(x) * 100)))
}

# Convert to long format
bar_data <- do.call(rbind, lapply(1:length(categories), function(i) {
  data.frame(
    freq = p_data[[categories[i]]],
    sample = rownames(p_data),
    category = categories[i]
  )
}))

# Dendrogram
dend_data <- dendro_data(
  as.dendrogram(hclust(dist(p_data, method = "euclidean"), method = "complete")),
  type = "rectangle"
)

segment_data <- dend_data$segments
label_data <- dend_data$labels

# Merge bar data with dendrogram positions
bar_data <- merge(bar_data, label_data, by.x = "sample", by.y = "label")

n_samples <- nrow(p_data)

# Layout parameters
dend_width <- 0.08
name_width <- 0.12
bar_len <- max(rowSums(p_data))
total_len <- bar_len / (1 - dend_width - name_width)

dend_max <- max(segment_data$y)
dend_len <- total_len * dend_width
segment_data$y <- segment_data$y / dend_max * dend_len
segment_data$yend <- segment_data$yend / dend_max * dend_len

name_len <- total_len * name_width

# Reference lines
y_breaks <- round(unlist(lapply(1:length(categories), function(i) {
  mean(rowSums(p_data[, 1:i, drop = FALSE]))
})), 1)

# Plot
p <- ggplot() +
  # Stacked bars
  geom_bar(
    data = bar_data,
    aes(x = x, y = freq, fill = factor(category, levels = rev(categories))),
    stat = "identity", width = 0.8
  ) +
  scale_fill_manual(values = cat_colors) +
  
  # Dendrogram
  geom_segment(
    data = segment_data,
    aes(x = x, y = -y - name_len, xend = xend, yend = -yend - name_len)
  ) +
  
  # Reference lines
  geom_segment(
    aes(x = 0.5, xend = n_samples + 0.5, y = y_breaks, yend = y_breaks),
    linetype = "dashed", alpha = 0.5
  ) +
  
  # Sample names
  geom_text(
    data = label_data,
    aes(x = x, y = -name_len / 2, label = label),
    size = 2, fontface = "bold", angle = 0
  ) +
  
  # Break labels
  geom_text(
    aes(x = 0.5, y = y_breaks, label = y_breaks),
    size = 2.5, fontface = "bold", hjust = 1
  ) +
  
  coord_flip() +
  labs(fill = "Category", title = "Gene Family Composition by Sample") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.position = "right"
  )

fig_width <- max(10, n_samples * 0.25)
ggsave(paste0(out_prefix, ".stackbar.pdf"), p, width = fig_width, height = 8, dpi = 300)
ggsave(paste0(out_prefix, ".stackbar.png"), p, width = fig_width, height = 8, dpi = 300)
cat("Saved:", paste0(out_prefix, ".stackbar.pdf/png"), "\n")
