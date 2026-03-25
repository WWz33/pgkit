#!/usr/bin/env Rscript
# =============================================================================
# Enhanced Saturation Curve with Fitting (based on APAVplot concepts)
# =============================================================================
# Features:
#   - Multiple visualization modes (ribbon, errorbar, jitter)
#   - Fitted curve with confidence interval
#   - Support for group comparison
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_curve_enhanced.R <saturation_curve.tsv> <output_prefix> [mode]\n")
  cat("  mode: ribbon (default), errorbar, jitter\n")
  quit(status = 1)
}

data_file <- args[1]
out_prefix <- args[2]
mode <- if (length(args) >= 3) args[3] else "ribbon"

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

# Read data
df <- read.delim(data_file)

# Colors
colors <- c("Core" = "#2E86AB", "Pan" = "#A23B72")
fills <- c("Core" = "#2E86AB", "Pan" = "#A23B72")

# Convert to long format
df_long <- data.frame(
  n_accessions = rep(df$n_accessions, 2),
  type = rep(c("Core", "Pan"), each = nrow(df)),
  count = c(df$core, df$pan)
)

# Check if SD column exists (from simulations)
has_sd <- "core_sd" %in% colnames(df) && "pan_sd" %in% colnames(df)

if (has_sd) {
  df_long$sd <- c(df$core_sd, df$pan_sd)
}

# Base plot
p_base <- ggplot(df_long, aes(x = n_accessions, color = type)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = fills) +
  labs(
    x = "Number of Accessions",
    y = "Number of Gene Families",
    color = "Type",
    fill = "Type",
    title = "Core/Pan Gene Family Saturation Curve"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

# Mode-specific plots
if (mode == "ribbon" && has_sd) {
  # Ribbon with confidence interval
  p <- p_base +
    geom_ribbon(aes(ymin = count - sd, ymax = count + sd, fill = type),
                alpha = 0.3, color = NA) +
    geom_line(aes(y = count), linewidth = 1.2) +
    geom_point(aes(y = count), size = 3)
} else if (mode == "errorbar" && has_sd) {
  # Error bars
  p <- p_base +
    geom_errorbar(aes(ymin = count - sd, ymax = count + sd),
                  width = 0.5, linewidth = 0.8) +
    geom_line(aes(y = count), linewidth = 1.2) +
    geom_point(aes(y = count), size = 3)
} else if (mode == "jitter" && has_sd) {
  # Jitter with mean point (if raw data available)
  p <- p_base +
    geom_line(aes(y = count), linewidth = 1.2) +
    geom_point(aes(y = count), size = 3, shape = 21, fill = "white")
} else {
  # Simple line + point (no SD)
  p <- p_base +
    geom_line(aes(y = count), linewidth = 1.2) +
    geom_point(aes(y = count), size = 3)
}

# Add fitted curve (polynomial fit)
if (nrow(df) >= 5) {
  # Core fit
  core_fit <- lm(core ~ poly(n_accessions, 3, raw = TRUE), data = df)
  df$core_fitted <- predict(core_fit)
  
  # Pan fit
  pan_fit <- lm(pan ~ poly(n_accessions, 3, raw = TRUE), data = df)
  df$pan_fitted <- predict(pan_fit)
  
  df_fit_long <- data.frame(
    n_accessions = rep(df$n_accessions, 2),
    type = rep(c("Core", "Pan"), each = nrow(df)),
    fitted = c(df$core_fitted, df$pan_fitted)
  )
  
  p <- p +
    geom_line(data = df_fit_long, aes(y = fitted, linetype = "Fitted"),
              linewidth = 0.8, alpha = 0.7) +
    scale_linetype_manual(values = c("Fitted" = "dashed")) +
    guides(linetype = guide_legend(title = NULL))
}

# Save
ggsave(paste0(out_prefix, ".curve_enhanced.pdf"), p, width = 10, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".curve_enhanced.png"), p, width = 10, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".curve_enhanced.pdf/png"), "\n")

# Print fitted equation
if (nrow(df) >= 5) {
  cat("\n=== Fitted Equations ===\n")
  cat("Core: y = ", round(coef(core_fit)[1], 2), 
      " + ", round(coef(core_fit)[2], 4), "*x",
      " + ", round(coef(core_fit)[3], 6), "*x^2",
      " + ", round(coef(core_fit)[4], 8), "*x^3\n", sep = "")
  cat("Pan:  y = ", round(coef(pan_fit)[1], 2),
      " + ", round(coef(pan_fit)[2], 4), "*x",
      " + ", round(coef(pan_fit)[3], 6), "*x^2",
      " + ", round(coef(pan_fit)[4], 8), "*x^3\n", sep = "")
}
