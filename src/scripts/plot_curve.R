#!/usr/bin/env Rscript
# =============================================================================
# Pan/Core Saturation Curve with Fitting
# =============================================================================
# Models:
#   Pan-genome (Heaps' law):  pan = P1 * n^gamma + P2
#   Core-genome (exponential): core = C1 * exp(-C2 * n) + C3
#
# Usage: Rscript plot_curve.R <saturation_curve.tsv> <output_prefix>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_curve.R <saturation_curve.tsv> <output_prefix>\n")
  quit(status = 1)
}

data_file <- args[1]
out_prefix <- args[2]

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# Read data
df <- read.delim(data_file)

n_max <- max(df$n_accessions)
pan_max <- max(df$pan)
core_final <- df$core[nrow(df)]

# ============================================================
# Heaps' law fitting for Pan-genome: pan = P1 * n^gamma + P2
# ============================================================
heaps_law <- function(n, P1, gamma, P2) P1 * n^gamma + P2

tryCatch({
  fit_pan <- nls(pan ~ heaps_law(n_accessions, P1, gamma, P2),
                 data = df,
                 start = list(P1 = pan_max * 0.5, gamma = 0.5, P2 = pan_max * 0.5),
                 control = nls.control(maxiter = 200))
  pan_coef <- coef(fit_pan)
  pan_fitted <- TRUE
  pan_eq <- sprintf("Pan = %.1f * n^%.3f + %.1f", pan_coef["P1"], pan_coef["gamma"], pan_coef["P2"])
  cat("Pan-genome fit:", pan_eq, "\n")
}, error = function(e) {
  # Fallback: simple power fit
  fit_pan_log <- lm(log(pan) ~ log(n_accessions), data = df[-1, ])
  pan_fitted <<- TRUE
  pan_coef <<- c(P1 = exp(coef(fit_pan_log)[1]), gamma = coef(fit_pan_log)[2], P2 = 0)
  pan_eq <<- sprintf("Pan = %.1f * n^%.3f", pan_coef["P1"], pan_coef["gamma"])
  cat("Pan-genome fit (log-linear):", pan_eq, "\n")
})

# ============================================================
# Exponential decay fitting for Core-genome: core = C1 * exp(-C2 * n) + C3
# ============================================================
exp_decay <- function(n, C1, C2, C3) C1 * exp(-C2 * n) + C3

tryCatch({
  fit_core <- nls(core ~ exp_decay(n_accessions, C1, C2, C3),
                  data = df,
                  start = list(C1 = core_final, C2 = 0.01, C3 = core_final * 0.8),
                  control = nls.control(maxiter = 200))
  core_coef <- coef(fit_core)
  core_fitted <- TRUE
  core_eq <- sprintf("Core = %.1f * exp(-%.4f * n) + %.1f", core_coef["C1"], core_coef["C2"], core_coef["C3"])
  cat("Core-genome fit:", core_eq, "\n")
}, error = function(e) {
  # Fallback: linear decay
  fit_core_log <- lm(log(max(core) - core) ~ n_accessions, data = df[-1, ])
  core_fitted <<- TRUE
  core_coef <<- c(C1 = core_final, C2 = -coef(fit_core_log)[2], C3 = 0)
  core_eq <<- sprintf("Core = %.1f (linear)", core_final)
  cat("Core-genome fit (linear):", core_eq, "\n")
})

# Generate fitted curves
n_seq <- seq(1, n_max, length.out = 100)

if (exists("fit_pan")) {
  pan_curve <- data.frame(n_accessions = n_seq,
                          count = predict(fit_pan, newdata = data.frame(n_accessions = n_seq)),
                          type = "Pan (fitted)")
} else {
  pan_curve <- data.frame(n_accessions = n_seq,
                          count = heaps_law(n_seq, pan_coef["P1"], pan_coef["gamma"], pan_coef["P2"]),
                          type = "Pan (fitted)")
}

if (exists("fit_core")) {
  core_curve <- data.frame(n_accessions = n_seq,
                           count = predict(fit_core, newdata = data.frame(n_accessions = n_seq)),
                           type = "Core (fitted)")
} else {
  core_curve <- data.frame(n_accessions = n_seq,
                           count = exp_decay(n_seq, core_coef["C1"], core_coef["C2"], core_coef["C3"]),
                           type = "Core (fitted)")
}

# Original data
df_plot <- data.frame(
  n_accessions = rep(df$n_accessions, 2),
  type = rep(c("Core", "Pan"), each = nrow(df)),
  count = c(df$core, df$pan)
)

# Confidence ribbon (if SD available)
has_sd <- "core_sd" %in% colnames(df) && "pan_sd" %in% colnames(df)

# Colors
colors <- c("Core" = "#2E86AB", "Pan" = "#A23B72",
            "Core (fitted)" = "#2E86AB", "Pan (fitted)" = "#A23B72")

# Plot
p <- ggplot()

# Add ribbon if SD available
if (has_sd) {
  df_ribbon <- data.frame(
    n_accessions = rep(df$n_accessions, 2),
    type = rep(c("Core", "Pan"), each = nrow(df)),
    count = c(df$core, df$pan),
    ymin = c(df$core - df$core_sd, df$pan - df$pan_sd),
    ymax = c(df$core + df$core_sd, df$pan + df$pan_sd)
  )
  p <- p +
    geom_ribbon(data = df_ribbon,
                aes(x = n_accessions, ymin = ymin, ymax = ymax, fill = type),
                alpha = 0.2, color = NA) +
    scale_fill_manual(values = colors)
}

# Add original lines (no fitted curves)
colors <- c("Core" = "#2E86AB", "Pan" = "#A23B72")

p <- p +
  geom_line(data = df_plot, aes(x = n_accessions, y = count, color = type),
            linewidth = 1.2) +

  # Scale and theme
  scale_color_manual(values = colors) +
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
    legend.position = "right"
  )

# Save
ggsave(paste0(out_prefix, ".curve.pdf"), p, width = 10, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".curve.png"), p, width = 10, height = 6, dpi = 300)
cat("Saved:", paste0(out_prefix, ".curve.pdf/png"), "\n")

# Output equations
cat("\n=== Fitted Equations ===\n")
cat(pan_eq, "\n")
cat(core_eq, "\n")
