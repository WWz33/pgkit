#!/usr/bin/env Rscript
# =============================================================================
# pgkit Color Palette Module
# =============================================================================
# Usage: source("palette.R")
# =============================================================================

# ============================================================
# Default color scheme
# ============================================================
PGKIT_PALETTE <- c(
  "#f78d85",  # 1 - coral pink
  "#ffd7cf",  # 2 - light pink
  "#48b6a6",  # 3 - teal
  "#b6dfd7",  # 4 - light teal
  "#85a6af",  # 5 - steel blue
  "#8d9dc7",  # 6 - lavender blue
  "#ffc725",  # 7 - gold
  "#ef7c11"   # 8 - orange
)

# Additional colors for extended palette
PGKIT_PALETTE_EXTENDED <- c(
  PGKIT_PALETTE,
  "#e74c3c",  # 9 - red
  "#3498db",  # 10 - blue
  "#2ecc71",  # 11 - green
  "#9b59b6",  # 12 - purple
  "#1abc9c",  # 13 - turquoise
  "#e67e22",  # 14 - dark orange
  "#34495e",  # 15 - dark gray
  "#f39c12"   # 16 - amber
)

# ============================================================
# Category colors (selected from palette)
# ============================================================
# Core: coral pink (#f78d85)
# Softcore: gold (#ffc725)
# Dispensable: teal (#48b6a6)
# Private: lavender blue (#8d9dc7)

CATEGORY_COLORS <- c(
  Core       = "#f78d85",
  Softcore   = "#ffc725",
  Dispensable = "#48b6a6",
  Private   = "#8d9dc7"
)

# Alternative (lowercase keys)
category_colors <- c(
  core       = "#f78d85",
  softcore   = "#ffc725",
  dispensable = "#48b6a6",
  specific   = "#8d9dc7"
)

# ============================================================
# PAV colors
# ============================================================
PAV_COLORS <- c(
  absence  = "#F0F0F0",
  presence = "#48b6a6"
)

# ============================================================
# Functions
# ============================================================

#' Get category colors
#' @param names Category names
#' @return Named vector of colors
get_category_colors <- function(names = NULL) {
  if (is.null(names)) {
    return(CATEGORY_COLORS)
  }
  # Handle both uppercase and lowercase
  names_lower <- tolower(names)
  result <- category_colors[names_lower]
  names(result) <- names
  return(result)
}

#' Get colors from palette
#' @param n Number of colors
#' @param palette Color palette (default: pgkit)
#' @return Vector of colors
get_palette <- function(n, palette = "pgkit") {
  if (palette == "pgkit") {
    colors <- PGKIT_PALETTE
  } else if (palette == "extended") {
    colors <- PGKIT_PALETTE_EXTENDED
  } else {
    colors <- PGKIT_PALETTE
  }
  
  if (n <= length(colors)) {
    return(colors[1:n])
  } else {
    # Recycle colors
    return(rep(colors, length.out = n))
  }
}

#' Get PAV colors
#' @return Named vector of colors
get_pav_colors <- function() {
  return(PAV_COLORS)
}

# ============================================================
# Print palette info
# ============================================================
print_palette <- function() {
  cat("pgkit Color Palette\n")
  cat("===================\n\n")
  
  cat("Default palette:\n")
  for (i in seq_along(PGKIT_PALETTE)) {
    cat(sprintf("  %d: %s\n", i, PGKIT_PALETTE[i]))
  }
  
  cat("\nCategory colors:\n")
  for (cat in names(CATEGORY_COLORS)) {
    cat(sprintf("  %s: %s\n", cat, CATEGORY_COLORS[cat]))
  }
  
  cat("\nPAV colors:\n")
  for (pav in names(PAV_COLORS)) {
    cat(sprintf("  %s: %s\n", pav, PAV_COLORS[pav]))
  }
}

# Auto-print when sourced
# print_palette()
