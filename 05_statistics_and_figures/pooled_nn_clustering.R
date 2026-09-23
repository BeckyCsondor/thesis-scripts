## ------------------------------------------------------------------
## Pooled nearest-neighbor clustering test: do hexagonal nucleation
## points sit spatially closer together than a random subset of
## fibrils of the same size, drawn from the same local density?
##
## Each of the 34 real wheel-level nearest-neighbor distances is
## measured only against OTHER wheels within its own view (tomogram /
## orientation group); the 6800 null draws are built the same way, per
## view, then pooled together - raw coordinates from different views
## are never mixed, only the resulting distance measurements are.
##
## Mann-Whitney U (Python, scipy): p=1.72e-06, R (Clark-Evans-style,
## observed/null mean) = 0.65. Significance shown here via a bracket +
## asterisks rather than the raw p-value.
##
## Input: pooled_nn_data.csv - columns: group, value
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

RED       <- "#B96965"
RED_DARK  <- "#7A3D3A"
GREY      <- "#cccccc"
GREY_DARK <- "#888888"

df <- read_csv("pooled_nn_data.csv", show_col_types = FALSE)
df$group <- factor(df$group, levels = c("Observed", "Null"))

## SEM (not SD) - matches what's plotted: error bars show how precisely
## the MEAN is known, not the full spread of individual values. This
## makes the null bar's error bar look small since n=6800, which is
## expected and correct for SEM.
stats_df <- df %>%
  group_by(group) %>%
  summarise(mean_val = mean(value), sem = sd(value) / sqrt(n()), n = n())

## Recompute the significance test directly in R too (Mann-Whitney U),
## so the script is self-contained rather than relying on a hardcoded
## p-value from the earlier Python analysis
observed_vals <- df$value[df$group == "Observed"]
null_vals <- df$value[df$group == "Null"]
wtest <- wilcox.test(observed_vals, null_vals, alternative = "less")

sig_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}
stars <- sig_stars(wtest$p.value)

bracket_y <- max(stats_df$mean_val + stats_df$sem) + 10
tick_h <- 3

build_plot <- function(font_family) {
  ggplot(stats_df, aes(x = group, y = mean_val, fill = group)) +
    geom_col(color = c(RED_DARK, GREY_DARK), width = 0.6, linewidth = 1.1) +
    geom_errorbar(aes(ymin = mean_val - sem, ymax = mean_val + sem),
                  width = 0.18, linewidth = 1.0, color = "black") +
    scale_fill_manual(values = c(Observed = RED, Null = GREY), guide = "none") +
    ## significance bracket with stars, spanning both bars
    annotate("segment", x = 1, xend = 1, y = bracket_y - tick_h, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("segment", x = 2, xend = 2, y = bracket_y - tick_h, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("text", x = 1.5, y = bracket_y + 3, label = stars,
             fontface = "bold", size = 7, family = font_family) +
    scale_x_discrete(labels = c(Observed = "Observed\n(real wheels)",
                                  Null = "Null\n(random fibril\nsubsets)")) +
    scale_y_continuous(limits = c(0, bracket_y + 12), expand = c(0, 0)) +
    labs(
      title = sprintf("Pooled across all views\n(n=%d wheels, n=%d null draws)",
                       stats_df$n[stats_df$group == "Observed"],
                       stats_df$n[stats_df$group == "Null"]),
      x = NULL,
      y = "nearest-neighbor distance (nm)"
    ) +
    theme_minimal(base_size = 13, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, margin = margin(b = 10)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(4, "pt"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10.5, color = "black"),
      axis.title.y = element_text(size = 11, margin = margin(r = 8))
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
  safe_ggsave("pooled_nn_clustering.png", p, width = 4.6, height = 6.6, dpi = 300,
              bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  safe_ggsave("pooled_nn_clustering.png", p, width = 4.6, height = 6.6, dpi = 300, bg = "white")
}

safe_ggsave("pooled_nn_clustering.pdf", p, width = 4.6, height = 6.6, bg = "white")

message(sprintf("\nMann-Whitney U p-value = %.3e (%s)", wtest$p.value, stars))
