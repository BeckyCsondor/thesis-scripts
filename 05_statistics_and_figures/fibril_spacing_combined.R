## ------------------------------------------------------------------
## Fibril spacing: minimum approach vs. sustained parallel gap
##
## Three violins:
##   1. Minimum Distance (Pooled) - line-to-line nearest-neighbor gap,
##      minus 6 nm fibril width. 7 fibril pairs visually confirmed in
##      IMOD as fragmented single-fibril traces (not genuine close
##      contacts) excluded; upper-IQR outliers also excluded. Pos022
##      excluded throughout (tracing not considered reliable).
##   2. Sustained Gap (Pooled, Parallel) - restricted to fibril pairs
##      within 8 degrees of parallel and >=50% length overlap; gap
##      averaged across the full overlapping run (not just the closest
##      point). Nearest-neighbor search unrestricted by orientation
##      group.
##   3. Sustained Gap (Orientation Restricted) - same sustained-gap
##      method, but the nearest-neighbor search is restricted to
##      fibrils sharing the same orientation group.
##
## The "ns" bracket compares the two SUSTAINED GAP groups only (not the
## minimum-distance group, which measures a different, non-comparable
## quantity) - Mann-Whitney U, matching the Python analysis exactly.
##
## Input: spacing_three_way.csv - columns: group, gap_nm
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

MIN_COL <- "#6b8ba4"        # blue-grey for minimum distance
SUSTAINED_COL <- "#B96965"  # muted red for sustained gap (tau fibrils)

df <- read_csv("spacing_three_way.csv", show_col_types = FALSE)

group_order <- c("Minimum Distance (Pooled)",
                  "Sustained Gap (Pooled, Parallel)",
                  "Sustained Gap (Orientation Restricted)")
df$group <- factor(df$group, levels = group_order)

## separate 2-level variable just for fill color/legend, since the violins
## use a 3-level x-position but only 2 distinct colors (the two "sustained
## gap" groups share a color) - mapping fill directly to `group` would
## give a redundant 3-entry legend
df$metric <- factor(ifelse(grepl("^Minimum", df$group), "Minimum distance (closest point)",
                            "Sustained gap (parallel pairs, averaged over overlap)"),
                     levels = c("Minimum distance (closest point)",
                                "Sustained gap (parallel pairs, averaged over overlap)"))

stats_df <- df %>%
  group_by(group) %>%
  summarise(mean_gap = mean(gap_nm), median_gap = median(gap_nm),
            max_gap = max(gap_nm), n = n())

## Mann-Whitney U comparing ONLY the two sustained-gap groups (positions
## 2 and 3) - the minimum-distance group measures a different quantity
## and isn't a meaningful comparison target for this test
sustained_pooled_vals <- df$gap_nm[df$group == "Sustained Gap (Pooled, Parallel)"]
sustained_within_vals <- df$gap_nm[df$group == "Sustained Gap (Orientation Restricted)"]
mw <- wilcox.test(sustained_pooled_vals, sustained_within_vals, alternative = "two.sided")

## geometry for the significance bracket, spanning x=2 to x=3 only
bracket_y <- max(stats_df$max_gap) + 5.5
tick_h <- 0.6

build_plot <- function(font_family) {
  ggplot(df, aes(x = group, y = gap_nm, fill = metric)) +
    geom_violin(color = "#333333", alpha = 0.75, linewidth = 1.0,
                width = 0.6, trim = TRUE) +
    scale_fill_manual(values = setNames(c(MIN_COL, SUSTAINED_COL), levels(df$metric))) +
    ## median: solid black line, scoped to each violin's own width
    stat_summary(fun = median, geom = "crossbar", width = 0.5,
                 fatten = 1, color = "#1a1a1a", linewidth = 0.9) +
    ## mean: dotted line, same width, drawn per group from stats_df
    geom_segment(data = stats_df,
                 aes(x = as.numeric(group) - 0.25, xend = as.numeric(group) + 0.25,
                     y = mean_gap, yend = mean_gap),
                 inherit.aes = FALSE, color = "#333333", linewidth = 0.9,
                 linetype = "dotted") +
    geom_text(data = stats_df,
              aes(x = group, y = max_gap + 1.2, label = sprintf("Median: %.1f", median_gap)),
              inherit.aes = FALSE, fontface = "bold.italic", size = 3.4,
              color = "#1a1a1a", family = font_family) +
    ## significance bracket between the two SUSTAINED GAP groups (x=2, x=3)
    annotate("segment", x = 2, xend = 2, y = bracket_y - tick_h, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("segment", x = 3, xend = 3, y = bracket_y - tick_h, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("segment", x = 2, xend = 3, y = bracket_y, yend = bracket_y,
             linewidth = 0.9, color = "black") +
    annotate("text", x = 2.5, y = bracket_y + 0.6, label = "ns",
             fontface = "bold", size = 4.2, family = font_family) +
    scale_y_continuous(limits = c(0, bracket_y + 4)) +
    scale_x_discrete(labels = sprintf("%s\n(n=%d)", gsub(" \\(", "\n(", stats_df$group), stats_df$n)) +
    labs(
      title = "Fibril spacing: minimum approach vs. sustained parallel gap",
      x = NULL,
      y = "Fibril Surface-Surface Distance (nm)",
      caption = paste0(
        "Minimum distance: line-to-line nearest-neighbor gap minus 6 nm fibril width; 7 pairs visually\n",
        "confirmed in IMOD as fragmented single-fibril traces excluded. Sustained gap: averaged across\n",
        "the full overlap of pairs within 8\u00b0 of parallel and \u226550% length overlap. Pos022 excluded throughout.\n",
        sprintf("ns: Mann-Whitney U (sustained gap groups only), p=%.3f", mw$p.value)
      )
    ) +
    theme_minimal(base_size = 12, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(4, "pt"),
      axis.text.x = element_text(size = 9.3, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, margin = margin(r = 8)),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.3),
      legend.key.size = unit(10, "pt"),
      legend.margin = margin(b = 6),
      plot.caption = element_text(size = 7.6, color = "#444444", face = "italic",
                                   hjust = 0.5, margin = margin(t = 14)),
      plot.caption.position = "plot"
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
  safe_ggsave("fibril_spacing_combined.png", p, width = 7.2, height = 6.6, dpi = 300,
              bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  safe_ggsave("fibril_spacing_combined.png", p, width = 7.2, height = 6.6, dpi = 300, bg = "white")
}

safe_ggsave("fibril_spacing_combined.pdf", p, width = 7.2, height = 6.6, bg = "white")
