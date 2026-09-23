## ------------------------------------------------------------------
## Fibril spacing (overall nearest-neighbor) - box-and-whisker
##
## Data: line-to-line (segment-to-segment) nearest-neighbor distance,
## minus 6 nm fibril width, with gaps <=6 nm (likely mis-picked) and
## upper-IQR outliers excluded. Pos022 excluded (tracing not considered
## reliable). See spacing_overall.csv for the underlying per-fibril values.
##
## Input: spacing_overall.csv - columns: tomogram, gap_nm
## ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(readr)

FONT_FAMILY <- "Helvetica Neue"  # change to "Helvetica" if unavailable
use_ragg <- requireNamespace("ragg", quietly = TRUE)

df <- read_csv("spacing_overall.csv", show_col_types = FALSE)

tomo_order <- c("Pos018", "Pos019", "Pos023")
df$tomogram <- factor(df$tomogram, levels = tomo_order)

## build a "pooled" pseudo-group containing all tomograms combined,
## appended as an extra category on the same axis
pooled_df <- df %>% mutate(tomogram = "pooled")
plot_df <- bind_rows(df, pooled_df) %>%
  mutate(tomogram = factor(tomogram, levels = c(tomo_order, "pooled")))

pooled_mean <- mean(pooled_df$gap_nm)

## per-group mean, for the text labels above each box
means_df <- plot_df %>%
  group_by(tomogram) %>%
  summarise(mean_gap = mean(gap_nm), max_gap = max(gap_nm), n = n())

TOMO_FILL   <- "#4a6fa5"
POOLED_FILL <- "#6a4c93"
fill_map <- setNames(c(rep(TOMO_FILL, length(tomo_order)), POOLED_FILL),
                      c(tomo_order, "pooled"))

kw <- kruskal.test(gap_nm ~ tomogram, data = df)  # test on the 3 tomograms only

p <- ggplot(plot_df, aes(x = tomogram, y = gap_nm, fill = tomogram)) +
  geom_boxplot(width = 0.55, linewidth = 0.7, color = "#2c4a70",
               outlier.shape = 21, outlier.size = 1.8, outlier.color = "#555555",
               outlier.fill = "white", outlier.alpha = 0.7) +
  scale_fill_manual(values = fill_map, guide = "none") +
  geom_hline(yintercept = pooled_mean, color = "#c0392b", linewidth = 0.9,
             linetype = "dashed") +
  geom_text(data = means_df, aes(x = tomogram, y = max_gap + 1.2,
                                  label = sprintf("%.1f nm", mean_gap)),
            inherit.aes = FALSE, fontface = "bold", size = 3.6,
            color = "#333333", family = FONT_FAMILY) +
  annotate("text", x = 1, y = pooled_mean - 1.5,
           label = sprintf("Pooled mean = %.1f nm", pooled_mean),
           hjust = 0, size = 3.2, color = "#c0392b", family = FONT_FAMILY) +
  scale_x_discrete(labels = setNames(sprintf("%s\n(n=%d)", means_df$tomogram, means_df$n),
                                      as.character(means_df$tomogram))) +
  labs(
    title = "Fibril spacing (Pos022 excluded, 6 nm width correction)",
    x = NULL,
    y = "surface-to-surface gap (nm)",
    caption = sprintf(
      "Kruskal\u2013Wallis across the 3 tomograms: H=%.1f, p=%.4f\nLine-to-line nearest-neighbor distance minus 6 nm fibril width; gaps \u22646 nm\n(likely mis-picked) and IQR-outliers excluded. Pos022 excluded.",
      kw$statistic, kw$p.value)
  ) +
  theme_minimal(base_size = 12, base_family = FONT_FAMILY) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(4, "pt"),
    axis.text = element_text(size = 10.5, color = "black"),
    axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
    plot.caption = element_text(size = 8.3, color = "#444444", face = "italic",
                                 hjust = 0.5, margin = margin(t = 14)),
    plot.caption.position = "plot"
  )

print(p)

if (use_ragg) {
  ggsave("fibril_spacing_overall.png", p, width = 6.5, height = 6.0, dpi = 300,
         bg = "white", device = ragg::agg_png)
} else {
  message("Package 'ragg' not found - install.packages('ragg') for better font rendering.")
  ggsave("fibril_spacing_overall.png", p, width = 6.5, height = 6.0, dpi = 300, bg = "white")
}

## base pdf() avoids the XQuartz/cairo dependency (see prior script notes);
## substitutes standard Helvetica for "Helvetica Neue" visually.
ggsave("fibril_spacing_overall.pdf", p, width = 6.5, height = 6.0, bg = "white")
