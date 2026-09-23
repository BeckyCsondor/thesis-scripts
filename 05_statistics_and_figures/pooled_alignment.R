## ------------------------------------------------------------------
## Fibril alignment: pooled S2 nematic order and Parallelism Index (PI)
##
## Both metrics are computed WITHIN each fibril's own orientation group
## (not against a single whole-tomogram director) - this avoids penalizing
## minority-group fibrils just for belonging to a second, differently-
## oriented population, which would otherwise look like "poor alignment"
## when it's really just a second real population.
##
## PI = mean(|cos(theta)|) between each fibril's axis and its own group's
## nematic director. S2 = mean(1.5*cos(theta)^2 - 0.5), the standard
## nematic order parameter (0 = random, 1 = perfect alignment).
##
## Colors deliberately kept in the red family (not blue) to visually
## group this figure with the "sustained distance" figures rather than
## the earlier blue "converging pair" figure.
##
## Input: alignment_pooled.csv - columns: tomogram, fibril_idx, S2, PI
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

RED_DARK  <- "#7A3D3A"
RED_MAIN  <- "#B96965"
RED_LIGHT <- "#d9a8a5"

df <- read_csv("alignment_pooled.csv", show_col_types = FALSE)

pooled <- df %>%
  summarise(S2 = mean(S2), PI = mean(PI), n = n()) %>%
  pivot_longer(cols = c(S2, PI), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("S2", "PI")))

n_fibrils <- nrow(df)

build_plot <- function(font_family) {
  ggplot(pooled, aes(x = metric, y = value, fill = metric)) +
    geom_col(width = 0.6, color = RED_DARK, linewidth = 1.0) +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.6,
              fontface = "bold", size = 4.2, family = font_family) +
    geom_hline(yintercept = 1.0, color = "gray50", linewidth = 0.6, linetype = "dotted") +
    scale_fill_manual(values = c(S2 = RED_MAIN, PI = RED_LIGHT),
                       labels = c(S2 = expression(S[2]~"order parameter"),
                                  PI = "Parallelism Index (PI)")) +
    scale_x_discrete(labels = c(S2 = expression(S[2]), PI = "PI")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +  # lets the value labels peek slightly
                                       # above y=1 without changing the axis
    labs(
      title = sprintf("Fibril alignment (within orientation group)\npooled, n=%d fibrils", n_fibrils),
      x = NULL,
      y = "Value (0 = random, 1 = perfect)"
    ) +
    theme_minimal(base_size = 13, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(4, "pt"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10),  # headroom for
                                                                # the labels
                                                                # peeking above y=1
      ## legend placed OUTSIDE the plotting area (right side) so it can
      ## never overlap the bars, regardless of bar height
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.size = unit(14, "pt")
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
  safe_ggsave("pooled_alignment.png", p, width = 4.6, height = 5.8, dpi = 300,
              bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  safe_ggsave("pooled_alignment.png", p, width = 4.6, height = 5.8, dpi = 300, bg = "white")
}

safe_ggsave("pooled_alignment.pdf", p, width = 4.6, height = 5.8, bg = "white")
