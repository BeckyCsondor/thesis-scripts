# ============================================================
#  PI/Hoechst Viability Analysis — R Plotting Script
#  Mouse hippocampal organotypic slices
#  Conditions: Control, wtTau, 2xTau
#
#  Output:
#  Figure A — dot plot with individual values and mean line
#  Figure B — bar graph of mean % dead nuclei per condition
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)   # for combining panels; install if needed

# ============================================================
#  DATA — segmentation-based % dead nuclei (PI+ / Hoechst)
#  From PI_Hoechst_results CSV, Seg_Pct_Dead column
# ============================================================

data <- data.frame(
  Sample    = c("Control_S1", "Control_S2",
                "wtTau_S1",   "wtTau_S2",
                "2xTau_S1",   "2xTau_S2"),
  Condition = factor(
    c("Control", "Control", "wtTau", "wtTau", "2xTau", "2xTau"),
    levels = c("Control", "wtTau", "2xTau")
  ),
  Pct_Dead  = c(1.6, 6.8, 7.3, 13.0, 36.1, 33.5)
)

# Condition means for bar graph and mean line on dot plot
means <- data %>%
  group_by(Condition) %>%
  summarise(Mean_Dead = mean(Pct_Dead), .groups = "drop")

# Condition colours — consistent across both panels
cond_cols <- c(
  "Control" = "#6aab7f",
  "wtTau"   = "#d4893a",
  "2xTau"   = "#c0504d"
)

# ============================================================
#  FIGURE A — Dot plot with individual values and mean line
# ============================================================

fig_a <- ggplot() +
  # Individual data points
  geom_point(
    data     = data,
    aes(x    = Condition, y = Pct_Dead, colour = Condition),
    size     = 4,
    alpha    = 0.9
  ) +
  # Mean line per condition
  geom_errorbar(
    data     = means,
    aes(x    = Condition,
        ymin = Mean_Dead,
        ymax = Mean_Dead,
        colour = Condition),
    width    = 0.35,
    linewidth = 1.2
  ) +
  scale_colour_manual(values = cond_cols) +
  scale_y_continuous(
    limits = c(0, 45),
    breaks = seq(0, 45, by = 10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x       = NULL,
    y       = "% dead nuclei (PI+)",
    caption = "Each point = one organotypic slice.\nHorizontal line = condition mean. n = 2 per condition."
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "none",
    axis.line          = element_line(colour = "gray30"),
    axis.ticks         = element_line(colour = "gray30"),
    plot.caption       = element_text(size = 8, colour = "gray50", hjust = 0),
    panel.grid.major.y = element_line(colour = "gray90", linewidth = 0.3)
  )

# ============================================================
#  FIGURE B — Bar graph of mean % dead per condition
# ============================================================

fig_b <- ggplot(means, aes(x = Condition, y = Mean_Dead, fill = Condition)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = cond_cols) +
  scale_y_continuous(
    limits = c(0, 45),
    breaks = seq(0, 45, by = 10),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x       = NULL,
    y       = "Mean % dead nuclei (PI+)",
    caption = "Bar = condition mean. n = 2 per condition."
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "none",
    axis.line          = element_line(colour = "gray30"),
    axis.ticks         = element_line(colour = "gray30"),
    plot.caption       = element_text(size = 8, colour = "gray50", hjust = 0),
    panel.grid.major.y = element_line(colour = "gray90", linewidth = 0.3)
  )

# ============================================================
#  COMBINE AND SAVE
# ============================================================

# Print individually
print(fig_a)
print(fig_b)

# Combined panel using patchwork
combined <- fig_a + fig_b +
  plot_annotation(
    title   = "PI/Hoechst viability — mouse hippocampal organotypics",
    caption = "Segmentation-based quantification. Threshold = 15, size filter 10-500 px².",
    tag_levels = "A",
    theme = theme(
      plot.title   = element_text(size = 13, face = "plain"),
      plot.caption = element_text(size = 8, colour = "gray50")
    )
  )

print(combined)

# Save individual panels
ggsave("viability_dotplot.pdf",  plot = fig_a,    width = 4,   height = 5)
ggsave("viability_dotplot.png",  plot = fig_a,    width = 4,   height = 5, dpi = 300)
ggsave("viability_barplot.pdf",  plot = fig_b,    width = 4,   height = 5)
ggsave("viability_barplot.png",  plot = fig_b,    width = 4,   height = 5, dpi = 300)
ggsave("viability_combined.pdf", plot = combined, width = 8,   height = 5)
ggsave("viability_combined.png", plot = combined, width = 8,   height = 5, dpi = 300)

message("Saved: viability_dotplot, viability_barplot, viability_combined — PDF and PNG")
