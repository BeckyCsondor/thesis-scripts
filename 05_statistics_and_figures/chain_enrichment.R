## ------------------------------------------------------------------
## Chain membership enrichment at hexagonal nucleation points
##
## Three groups compared: baseline (all fibrils), automated hexagonal
## candidates (chain-blind geometric criterion, gap std < 20deg), and
## manually-annotated nucleation points. Colors deliberately chosen as
## a fresh pastel set that doesn't match any earlier figure in this
## thesis (avoids the tau-red B96965, the blue-grey convergence color,
## and the gold/orange nucleation-point color already used elsewhere).
##
## Input: validation_summary.csv - columns: group, n, pct_in_chain, pvalue_label
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

## professional pastel palette - muted sage, dusty lavender, soft apricot
PASTEL_SAGE    <- "#9CB89C"
PASTEL_LAVENDER <- "#B3A7CC"
PASTEL_APRICOT  <- "#E8B98A"
PASTEL_SAGE_DARK    <- "#5C7A5C"
PASTEL_LAVENDER_DARK <- "#6B5C8A"
PASTEL_APRICOT_DARK  <- "#A6702F"

df <- read_csv("validation_summary.csv", show_col_types = FALSE)
df$group <- factor(df$group, levels = c("All fibrils (baseline)",
                                          "Automated hexagonal candidates",
                                          "Manual nucleation points"))

fill_vals  <- c(PASTEL_SAGE, PASTEL_LAVENDER, PASTEL_APRICOT)
edge_vals  <- c(PASTEL_SAGE_DARK, PASTEL_LAVENDER_DARK, PASTEL_APRICOT_DARK)
names(fill_vals) <- levels(df$group)
names(edge_vals) <- levels(df$group)

build_plot <- function(font_family) {
  ggplot(df, aes(x = group, y = pct_in_chain, fill = group, color = group)) +
    geom_col(width = 0.6, linewidth = 1.1) +
    geom_text(aes(label = sprintf("%.0f%%", pct_in_chain)), vjust = -1.6,
              fontface = "bold", size = 5, family = font_family, color = "#333333") +
    geom_text(aes(label = pvalue_label), vjust = -3.4,
              fontface = "italic", size = 3.2, family = font_family, color = "#555555") +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_color_manual(values = edge_vals, guide = "none") +
    scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
    scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Two independent methods agree: hexagonal nucleation\npoints are strongly enriched for chain membership",
      x = NULL,
      y = "% that are members of a chain"
    ) +
    theme_minimal(base_size = 13, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 14.5, hjust = 0.5, margin = margin(b = 14)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(4, "pt"),
      axis.text.x = element_text(size = 10.5, color = "black"),
      axis.text.y = element_text(size = 10.5, color = "black"),
      axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
      plot.margin = margin(t = 25, r = 15, b = 10, l = 10)
    )
}

p <- build_plot(FONT_FAMILY)

tryCatch({
  print(p)
}, error = function(e) {
  message("On-screen preview failed (font-rendering quirk) - saved files below unaffected.")
})

safe_ggsave <- function(filename, plot_obj, ...) {
  tryCatch({
    ggsave(filename, plot_obj, ...)
    message("Saved ", filename, " using font '", FONT_FAMILY, "'.")
  }, error = function(e) {
    message("Saving with font '", FONT_FAMILY, "' failed - retrying with default font.")
    fallback_plot <- build_plot("")
    ggsave(filename, fallback_plot, ...)
    message("Saved ", filename, " using the default font as a fallback.")
  })
}

if (use_ragg) {
  safe_ggsave("chain_enrichment_pastel.png", p, width = 6.2, height = 6.6, dpi = 300,
              bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  safe_ggsave("chain_enrichment_pastel.png", p, width = 6.2, height = 6.6, dpi = 300, bg = "white")
}

safe_ggsave("chain_enrichment_pastel.pdf", p, width = 6.2, height = 6.6, bg = "white")
