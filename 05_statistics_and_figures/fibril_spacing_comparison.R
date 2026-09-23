## ------------------------------------------------------------------
## Fibril spacing: pooled (any neighbor) vs orientation-group-restricted
## (nearest neighbor within the same orientation group only) - violin
##
## Both distributions use the identical pipeline: line-to-line
## (segment-to-segment) nearest-neighbor distance, minus 6 nm fibril
## width, gaps <=10 nm excluded (visually confirmed in IMOD to be
## fragmented tracing of a single fibril, not genuine close packing),
## upper-IQR outliers excluded. Pos022 excluded throughout (tracing not
## considered reliable). Pos018 + Pos019 + Pos023 pooled in both cases.
##
## "Pooled Fibrils": nearest neighbor can be any fibril, any orientation
## group.
## "Orientation Group Restricted": nearest-neighbor search restricted to
## other fibrils in the SAME orientation group only (cross-group fibrils
## are never candidates, even if physically closer).
##
## Input: spacing_comparison.csv - columns: group, gap_nm
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

FIBRIL_COL <- "#B96965"  # muted red - represents tau fibrils

df <- read_csv("spacing_comparison.csv", show_col_types = FALSE)

group_order <- c("Pooled Fibrils", "Orientation Group Restricted")
df$group <- factor(df$group, levels = group_order)

stats_df <- df %>%
  group_by(group) %>%
  summarise(mean_gap = mean(gap_nm), median_gap = median(gap_nm),
            max_gap = max(gap_nm), n = n())

## Mann-Whitney U (location shift) and Kolmogorov-Smirnov (distribution
## shape) tests comparing the two groups, annotated on the plot
overall_vals <- df$gap_nm[df$group == "Pooled Fibrils"]
within_vals  <- df$gap_nm[df$group == "Orientation Group Restricted"]
mw <- wilcox.test(overall_vals, within_vals, alternative = "two.sided")
ks <- suppressWarnings(ks.test(overall_vals, within_vals))

## geometry for the significance bracket, placed above both violins
bracket_y <- max(stats_df$max_gap) + 5.5
tick_h <- 0.6

## Build the plot as a function of font family, so we can retry with a
## safe fallback if the requested font can't be rendered by the save
## device (this is what threw "invalid font type" previously).
build_plot <- function(font_family) {
  ggplot(df, aes(x = group, y = gap_nm)) +
    geom_violin(fill = FIBRIL_COL, color = "#7A3D3A", alpha = 0.75,
                linewidth = 1.0, width = 0.6, trim = TRUE) +
    ## median: solid black line, scoped to each violin's own width
    stat_summary(fun = median, geom = "crossbar", width = 0.5,
                 fatten = 1, color = "#1a1a1a", linewidth = 0.9) +
    ## mean: dotted line, same width, drawn per group from stats_df
    geom_segment(data = stats_df,
                 aes(x = as.numeric(group) - 0.25, xend = as.numeric(group) + 0.25,
                     y = mean_gap, yend = mean_gap),
                 inherit.aes = FALSE, color = "#3a1a18", linewidth = 0.9,
                 linetype = "dotted") +
    geom_text(data = stats_df,
              aes(x = group, y = max_gap + 1.2, label = sprintf("Median: %.1f", median_gap)),
              inherit.aes = FALSE, fontface = "bold.italic", size = 3.6,
              color = "#1a1a1a", family = font_family) +
    geom_text(data = stats_df,
              aes(x = as.numeric(group) + 0.32, y = mean_gap - 0.8,
                  label = sprintf("Mean: %.1f", mean_gap)),
              inherit.aes = FALSE, fontface = "italic", size = 3.3,
              color = "#7A3D3A", family = font_family, hjust = 0) +
    ## significance bracket (ns, since neither test reached p<0.05)
    annotate("segment", x = 1, xend = 1, y = bracket_y - tick_h, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("segment", x = 2, xend = 2, y = bracket_y - tick_h, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("text", x = 1.5, y = bracket_y + 0.6, label = "ns",
             fontface = "bold", size = 4.4, family = font_family) +
    scale_y_continuous(limits = c(0, bracket_y + 2.5)) +
    scale_x_discrete(labels = sprintf("%s\n(n=%d)", stats_df$group, stats_df$n)) +
    labs(
      title = "Surface-Surface Distance of Fibrils",
      x = NULL,
      y = "Fibril Surface-Surface Distance (nm)",
      caption = paste0(
        "Line-to-line nearest-neighbor distance minus 6 nm fibril width; gaps <=10 nm\n",
        "(visually confirmed fragmented traces) and IQR-outliers excluded. Pos022 excluded.\n",
        sprintf("Mann-Whitney U p=%.3f  |  Kolmogorov-Smirnov p=%.3f  -  no significant difference",
                mw$p.value, ks$p.value)
      )
    ) +
    theme_minimal(base_size = 12, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 15, margin = margin(b = 22)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(4, "pt"),
      axis.text = element_text(size = 10.5, color = "black"),
      axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
      plot.caption = element_text(size = 8.1, color = "#444444", face = "italic",
                                   hjust = 0.5, margin = margin(t = 14)),
      plot.caption.position = "plot"
    )
}

p <- build_plot(FONT_FAMILY)

## On-screen preview: the default interactive graphics device often
## can't resolve a named system font like "Helvetica Neue", throwing
## "invalid font type". This doesn't affect the saved files below.
tryCatch({
  print(p)
}, error = function(e) {
  message("On-screen preview failed (likely a font-rendering quirk of the ",
          "interactive graphics device) - this does NOT affect the saved ",
          "PNG/PDF files below. Error was: ", conditionMessage(e))
})

## Safe save: try with the requested font first; if ANY part of the
## rendering pipeline can't resolve it (this is what caused the
## previous crash, specifically on the "ns" bracket text), silently
## rebuild the plot with the system default font and save that instead,
## so a run always produces a usable file.
safe_ggsave <- function(filename, plot_obj, ...) {
  tryCatch({
    ggsave(filename, plot_obj, ...)
    message("Saved ", filename, " using font '", FONT_FAMILY, "'.")
  }, error = function(e) {
    message("Saving with font '", FONT_FAMILY, "' failed (", conditionMessage(e),
            ") - retrying with the default font instead.")
    fallback_plot <- build_plot("")  # "" = device default font
    ggsave(filename, fallback_plot, ...)
    message("Saved ", filename, " using the default font as a fallback.")
  })
}

if (use_ragg) {
  safe_ggsave("fibril_spacing_comparison.png", p, width = 5.2, height = 6.0, dpi = 300,
              bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  safe_ggsave("fibril_spacing_comparison.png", p, width = 5.2, height = 6.0, dpi = 300, bg = "white")
}

safe_ggsave("fibril_spacing_comparison.pdf", p, width = 5.2, height = 6.0, bg = "white")
