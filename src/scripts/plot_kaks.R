#!/usr/bin/env Rscript
# =============================================================================
# Ka/Ks Boxplot with Statistical Tests
# =============================================================================
# Methods:
#   1. Kruskal-Wallis test (non-parametric)
#   2. Fisher's LSD on ranks (alpha=0.001) - original paper method
#   3. Wilcoxon pairwise with BH correction - more robust
#
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
for (pkg in c("ggplot2", "agricolae", "ggpubr")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(ggplot2)
library(agricolae)
library(ggpubr)

# Default category colors
cat_colors <- c(
  Core = "#f78d85",
  Softcore = "#ffc725",
  Dispensable = "#48b6a6",
  Private = "#8d9dc7"
)

# Read data
df <- read.delim(kaks_file)
df <- df[!is.na(df$Ka_Ks) & is.finite(df$Ka_Ks) & df$Ka_Ks > 0 & df$Ka_Ks < 5, ]
df$Category <- factor(df$Category, levels = c("Core", "Softcore", "Dispensable", "Private"))
df <- df[!is.na(df$Category), ]

cat("=== Ka/Ks Analysis ===\n\n")

# Summary
cat("Summary:\n")
for (cat in levels(df$Category)) {
  sub <- df[df$Category == cat, ]
  if (nrow(sub) > 0) {
    cat(sprintf("  %s: n=%d, median=%.4f, mean=%.4f\n", cat, nrow(sub),
                median(sub$Ka_Ks), mean(sub$Ka_Ks)))
  }
}

# ============================================================
# Method 1: Kruskal-Wallis + Fisher's LSD (original paper)
# ============================================================
cat("\n--- Method 1: Fisher's LSD (alpha=0.001) ---\n")

kw <- kruskal.test(Ka_Ks ~ Category, data = df)
cat(sprintf("Kruskal-Wallis: chi-sq=%.4f, p=%.2e\n", kw$statistic, kw$p.value))

lsd_groups <- NULL
if (kw$p.value < 0.05) {
  df$rank_KaKs <- rank(df$Ka_Ks)
  aov_model <- aov(rank_KaKs ~ Category, data = df)
  lsd_result <- LSD.test(aov_model, "Category", p.adj = "none", alpha = 0.001)
  lsd_groups <- lsd_result$groups
  lsd_groups$Category <- rownames(lsd_groups)
  cat("\nFisher's LSD groups (alpha=0.001):\n")
  print(lsd_groups[, c("Category", "groups")])
  write.csv(lsd_groups, file = paste0(out_prefix, "_kaks_fisher_lsd.csv"), row.names = FALSE)
}

# ============================================================
# Method 2: Wilcoxon pairwise with BH correction (more robust)
# ============================================================
cat("\n--- Method 2: Wilcoxon pairwise (BH correction) ---\n")

pw <- pairwise.wilcox.test(df$Ka_Ks, df$Category, p.adjust.method = "BH")
cat("\nPairwise p-values:\n")
print(pw)

# Save pairwise results
pw_df <- as.data.frame(pw$p.value)
write.csv(pw_df, file = paste0(out_prefix, "_kaks_pairwise_wilcox.csv"))

# Get significance stars
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# ============================================================
# Boxplot with both annotations
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

# Add Fisher's LSD letters
if (!is.null(lsd_groups) && nrow(lsd_groups) > 0) {
  max_vals <- tapply(df$Ka_Ks, df$Category, max)
  annot_lsd <- data.frame(
    Category = lsd_groups$Category,
    y_pos = max_vals[lsd_groups$Category] * 1.08,
    label = lsd_groups$groups
  )
  annot_lsd <- annot_lsd[!is.na(annot_lsd$y_pos), ]
  
  if (nrow(annot_lsd) > 0) {
    p <- p + geom_text(
      data = annot_lsd,
      aes(x = Category, y = y_pos, label = label),
      size = 5, fontface = "bold"
    )
  }
}

# Add pairwise comparisons with stars
existing_cats <- levels(df$Category)
my_comparisons <- list(
  c("Core", "Softcore"),
  c("Core", "Dispensable"),
  c("Softcore", "Dispensable")
)
my_comparisons <- Filter(function(x) all(x %in% existing_cats), my_comparisons)

if (length(my_comparisons) > 0 && kw$p.value < 0.05) {
  p <- p + stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    tip.length = 0.02,
    step.increase = 0.08
  )
}

# Add global p-value
p <- p + annotate(
  "text",
  x = Inf, y = Inf,
  label = sprintf("Kruskal-Wallis p = %.2e", kw$p.value),
  hjust = 1.1, vjust = 1.5,
  size = 3.5, fontface = "italic"
)

# Save
ggsave(paste0(out_prefix, "_kaks_boxplot.pdf"), p, width = 8, height = 7, dpi = 300)
ggsave(paste0(out_prefix, "_kaks_boxplot.png"), p, width = 8, height = 7, dpi = 300)
cat("\nSaved:", paste0(out_prefix, "_kaks_boxplot.pdf/png"), "\n")

cat("\n=== Done ===\n")
