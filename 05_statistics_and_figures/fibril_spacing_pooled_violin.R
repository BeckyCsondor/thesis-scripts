## ------------------------------------------------------------------
## Fibril spacing (pooled) - violin plot, final version
##
## Data: line-to-line (segment-to-segment) nearest-neighbor distance,
## minus 6 nm fibril width. Gaps <=10 nm excluded - visual inspection
## in IMOD confirmed these represent fragmented tracing of a single
## fibril (should be merged), not genuinely distinct close-packed
## fibrils. Upper-IQR outliers also excluded. Pos022 excluded entirely
## (tracing not considered reliable). All remaining tomograms (Pos018,
## Pos019, Pos023) pooled into a single distribution.
##
## Input: spacing_overall.csv - columns: tomogram, gap_nm
## (already reflects the <=10nm exclusion - see accompanying notes)
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

FIBRIL_COL <- "#B96965"  # muted red - represents tau fibrils

df <- read_csv("spacing_overall.csv", show_col_types = FALSE)

pooled_df <- df %>% mutate(group = "pooled")
pooled_mean <- mean(pooled_df$gap_nm)
pooled_median <- median(pooled_df$gap_nm)
n_pooled <- nrow(pooled_df)
max_gap <- max(pooled_df$gap_nm)

p <- ggplot(pooled_df, aes(x = group, y = gap_nm)) +
  geom_violin(fill = FIBRIL_COL, color = "#7A3D3A", alpha = 0.75,
              linewidth = 1.0, width = 0.5, trim = TRUE) +
  ## median: solid black line, scoped to the violin width (stat_summary
  ## with a crossbar draws it only as wide as the violin, not the full axis)
  stat_summary(fun = median, geom = "crossbar", width = 0.42,
               fatten = 1, color = "#1a1a1a", linewidth = 0.9) +
  ## mean: dotted line, same width as the median crossbar
  annotate("segment", x = 1 - 0.21, xend = 1 + 0.21,
           y = pooled_mean, yend = pooled_mean,
           color = "#3a1a18", linewidth = 0.9, linetype = "dotted") +
  annotate("text", x = 1, y = max_gap + 1.2,
           label = sprintf("%.1f nm", pooled_mean),
           fontface = "bold", size = 4.2, color = "#333333", family = FONT_FAMILY) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_discrete(labels = sprintf("pooled\n(n=%d)", n_pooled),
                    expand = expansion(mult = 0.35)) +
  labs(
    title = "Fibril spacing\n(pooled, Pos022 excluded)",
    x = NULL,
    y = "surface-to-surface gap (nm)",
    caption = paste0(
      "Line-to-line nearest-neighbor distance minus 6 nm fibril width;\n",
      "gaps \u226410 nm (visually confirmed fragmented traces) and\n",
      "IQR-outliers excluded. Pos022 excluded."
    )
  ) +
  theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                               margin = margin(b = 10), lineheight = 1.05),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(4, "pt"),
    axis.text = element_text(size = 10.5, color = "black"),
    axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
    plot.caption = element_text(size = 8.1, color = "#444444", face = "italic",
                                 hjust = 0.5, margin = margin(t = 14)),
    plot.caption.position = "plot"
  )

print(p)

if (use_ragg) {
  ggsave("fibril_spacing_pooled_violin.png", p, width = 2.6, height = 5.8, dpi = 300,
         bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  ggsave("fibril_spacing_pooled_violin.png", p, width = 2.6, height = 5.8, dpi = 300, bg = "white")
}

ggsave("fibril_spacing_pooled_violin.pdf", p, width = 2.6, height = 5.8, bg = "white")
