#!/usr/bin/env Rscript
# =============================================================================
# Ka/Ks Boxplot with Kruskal-Wallis + Fisher's LSD Test
# =============================================================================
# Usage: Rscript plot_kaks.R <kaks_values.tsv> <output_prefix>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_kaks.R <kaks_values.tsv> <output_prefix>\n")
  quit(status = 1)
}

kaks_file <- args[1]
out_prefix <- args[2]

# Install packages if needed
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!require("agricolae", quietly = TRUE)) {
  install.packages("agricolae")
}
library(ggplot2)
library(agricolae)

# Default category colors
cat_colors <- c(
  Core = "#f78d85",
  Softcore = "#ffc725",
  Dispensable = "#48b6a6",
  Private = "#8d9dc7"
)

# Try to source palette.R from same directory as this script
script_path <- tryCatch({
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmdArgs)
  if (length(file_arg) > 0) {
    sub("--file=", "", cmdArgs[file_arg[1]])
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(script_path)) {
  palette_path <- file.path(dirname(script_path), "palette.R")
  if (file.exists(palette_path)) {
    source(palette_path)
  }
}

# Read data
df <- read.delim(kaks_file)

# Filter valid Ka/Ks values
df <- df[!is.na(df$Ka_Ks) & is.finite(df$Ka_Ks) & df$Ka_Ks > 0 & df$Ka_Ks < 5, ]

# Convert Category to factor
df$Category <- factor(df$Category, levels = c("Core", "Softcore", "Dispensable", "Private"))

# Remove empty categories
df <- df[!is.na(df$Category), ]

cat("=== Ka/Ks Analysis ===\n\n")

# ============================================================
# Summary statistics
# ============================================================
cat("Summary:\n")
for (cat in levels(df$Category)) {
  sub <- df[df$Category == cat, ]
  if (nrow(sub) > 0) {
    cat(sprintf("  %s: n=%d, median=%.4f, mean=%.4f\n", cat, nrow(sub),
                median(sub$Ka_Ks), mean(sub$Ka_Ks)))
  }
}

# ============================================================
# Kruskal-Wallis test
# ============================================================
annot <- NULL

if (length(unique(df$Category)) >= 2) {
  kw <- kruskal.test(Ka_Ks ~ Category, data = df)
  cat(sprintf("\nKruskal-Wallis test: chi-sq=%.4f, p=%.2e\n", kw$statistic, kw$p.value))
  
  # ============================================================
  # Fisher's LSD post-hoc test (alpha = 0.001)
  # Use kruskal + LSD.test with rank transformation
  # ============================================================
  if (kw$p.value < 0.05) {
    cat("\nPost-hoc test (Fisher's LSD, alpha=0.001):\n")
    
    # Rank transform for non-parametric post-hoc
    df$rank_KaKs <- rank(df$Ka_Ks)
    
    # Run LSD.test on ranks (need aov model)
    aov_model <- aov(rank_KaKs ~ Category, data = df)
    lsd_result <- LSD.test(aov_model, "Category", p.adj = "none", alpha = 0.001)
    
    # Get groups
    groups <- lsd_result$groups
    groups$Category <- rownames(groups)
    
    cat("\nGroups:\n")
    print(groups)
    
    # Save groups
    write.csv(groups, file = paste0(out_prefix, "_kaks_groups.csv"), row.names = FALSE)
    cat("\nGroups saved to:", paste0(out_prefix, "_kaks_groups.csv"), "\n")
    
    # Prepare annotation for boxplot
    max_vals <- tapply(df$Ka_Ks, df$Category, max)
    annot <- data.frame(
      Category = groups$Category,
      y_pos = max_vals[groups$Category] * 1.1,
      label = groups$groups
    )
    annot <- annot[!is.na(annot$y_pos), ]
  } else {
    cat("No significant difference (p >= 0.05), skipping post-hoc test\n")
  }
} else {
  cat("Only one category, skipping statistical tests\n")
}

# ============================================================
# Boxplot with significance letters
# ============================================================
p <- ggplot(df, aes(x = Category, y = Ka_Ks, fill = Category)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, outlier.alpha = 0.5) +
  scale_fill_manual(values = cat_colors) +
  labs(
    x = "Gene Category",
    y = "Ka/Ks",
    title = "Ka/Ks Distribution by Gene Category"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text = element_text(face = "bold")
  )

# Add significance letters if available
if (!is.null(annot) && nrow(annot) > 0) {
  p <- p + geom_text(
    data = annot,
    aes(x = Category, y = y_pos, label = label),
    size = 5, fontface = "bold", vjust = 0
  )
}

# Save
ggsave(paste0(out_prefix, "_kaks_boxplot.pdf"), p, width = 8, height = 6, dpi = 300)
ggsave(paste0(out_prefix, "_kaks_boxplot.png"), p, width = 8, height = 6, dpi = 300)
cat("\nSaved:", paste0(out_prefix, "_kaks_boxplot.pdf/png"), "\n")

cat("\n=== Done ===\n")
