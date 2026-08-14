# ============================================================
#  Nearest-Neighbour Distance Analysis
#  PI+ nuclei to Tau+ puncta — 2xTau hippocampal organotypics
#  
#  Output: histogram with distance distribution
# ============================================================

library(ggplot2)

# ============================================================
#  DATA — paste your distances from the Fiji Log here
# ============================================================

s2_distances <- c(
  5.47, 4.76, 6.87, 5.81, 0.49, 3.20, 8.48, 0.68, 4.65,
  10.23, 1.23, 3.68, 5.63, 1.27, 4.38, 5.83, 6.29, 0.23,
  2.72, 22.28, 14.68
)

s1_distances <- c(
  0.27, 4.07, 0.29, 48.61, 38.96, 16.41, 0.46, 0.43, 1.00,
  0.17, 0.49, 0.10, 0.38, 0.37, 0.35, 0.00, 0.37, 0.15,
  1.19, 0.31, 0.14, 13.61, 0.11, 5.65, 7.24, 0.29, 2.85,
  13.26, 5.17, 7.63, 1.38, 0.33, 0.12, 7.17, 13.04, 1.12,
  3.70, 2.55, 0.59, 1.88, 0.17, 0.33, 0.20, 0.32, 5.73,
  5.88, 0.36, 0.38, 0.06, 0.21, 11.11, 0.21, 0.76, 0.50,
  0.21, 0.13, 0.24, 4.67, 0.39, 1.75, 9.08, 0.16, 4.44,
  0.27, 0.14, 0.21, 1.10, 9.84, 0.09, 0.68, 0.29, 0.21,
  0.72, 0.62, 0.60, 1.53, 0.09, 0.34, 0.21, 0.26, 0.34,
  0.16, 3.50, 0.16, 0.20, 0.34, 0.32, 0.23, 5.02, 2.09,
  5.70, 0.21, 13.26, 0.20, 0.21, 0.09, 1.04, 0.24, 0.08,
  0.25, 14.95, 0.37, 0.20, 0.21, 1.46, 18.27, 0.40, 0.36,
  0.46, 4.26, 0.18, 3.06, 0.23, 0.25, 0.58, 9.61, 0.52,
  2.32, 0.23, 0.58, 0.55, 5.67, 0.61, 0.32, 4.44, 1.97,
  0.65, 0.50, 0.20, 0.09, 0.42, 0.30, 3.63, 0.18, 12.87,
  0.17, 3.23, 9.35, 15.37, 0.23, 0.27, 0.18, 17.66, 8.02,
  6.86, 2.43, 4.52, 2.46, 6.37, 8.94, 10.00, 0.24, 0.39,
  0.33, 3.39, 0.22
)

# Combine into one data frame
data <- data.frame(
  dist  = c(s1_distances, s2_distances),
  slice = c(rep("S1", length(s1_distances)), rep("S2", length(s2_distances)))
)

# Label whether within or beyond one cell diameter (~15 um)
cell_diameter <- 15

# Cap outliers for display
data$dist_capped <- pmin(data$dist, 50)

# Three proximity zones
same_cell_threshold     <- 5   # um - directly overlapping/same cell
adjacent_cell_threshold <- 15  # um - within one cell diameter

data$zone <- ifelse(data$dist <= same_cell_threshold,
                    "Same cell (0-5 \u00b5m)",
             ifelse(data$dist <= adjacent_cell_threshold,
                    "Adjacent cell (5-15 \u00b5m)",
                    "Not associated (>15 \u00b5m)"))

data$zone <- factor(data$zone, levels = c(
  "Same cell (0-5 \u00b5m)",
  "Adjacent cell (5-15 \u00b5m)",
  "Not associated (>15 \u00b5m)"
))

zone_cols <- c(
  "Same cell (0-5 \u00b5m)"        = "#c0504d",
  "Adjacent cell (5-15 \u00b5m)"   = "#d4897a",
  "Not associated (>15 \u00b5m)"   = "#e8c4bc"
)

# Summary stats
n_total    <- nrow(data)
n_same     <- sum(data$dist <= same_cell_threshold)
n_adjacent <- sum(data$dist > same_cell_threshold & data$dist <= adjacent_cell_threshold)
n_distal   <- sum(data$dist > adjacent_cell_threshold)
med_dist   <- round(median(data$dist), 2)

# ============================================================
#  PLOT
# ============================================================
p <- ggplot(data, aes(x = dist_capped, fill = zone)) +
  geom_histogram(
    binwidth  = 2.5,
    colour    = "white",
    linewidth = 0.3,
    boundary  = 0
  ) +
  geom_vline(
    xintercept = same_cell_threshold,
    linetype   = "dashed",
    colour     = "#c0504d",
    linewidth  = 0.6,
    alpha      = 0.7
  ) +
  geom_vline(
    xintercept = adjacent_cell_threshold,
    linetype   = "dashed",
    colour     = "#c0504d",
    linewidth  = 0.6,
    alpha      = 0.4
  ) +
  annotate("text",
           x = same_cell_threshold + 0.4, y = Inf,
           label = "same cell", hjust = 0, vjust = 1.8,
           size = 3, colour = "#c0504d") +
  annotate("text",
           x = adjacent_cell_threshold + 0.4, y = Inf,
           label = "~1 cell\ndiameter", hjust = 0, vjust = 1.3,
           size = 3, colour = "#c0504d", alpha = 0.7) +
  annotate("text",
           x = 35, y = Inf,
           label = paste0("Same cell: ", n_same, " (", round(n_same/n_total*100,1), "%)\n",
                          "Adjacent:  ", n_adjacent, " (", round(n_adjacent/n_total*100,1), "%)\n",
                          "Distal:    ", n_distal, " (", round(n_distal/n_total*100,1), "%)\n",
                          "Median = ", med_dist, " \u00b5m"),
           hjust = 0, vjust = 1.5, size = 3, colour = "gray30") +
  scale_fill_manual(values = zone_cols) +
  scale_x_continuous(
    limits = c(0, 55),
    breaks = seq(0, 50, by = 2.5),
    labels = function(x) {
      ifelse(x == 50, paste0(x, "+"), paste0(x, "-", x + 2.5))
    },
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x       = "Distance to nearest Tau+ punctum (\u00b5m)",
    y       = "Number of PI+ nuclei",
    fill    = NULL,
    caption = paste0(
      "2xTau hippocampal organotypic slices (S1 n=", length(s1_distances),
      ", S2 n=", length(s2_distances), "). ",
      "Bin width = 2.5 \u00b5m. Values >50 \u00b5m capped for display (n=2)."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 10),
    plot.caption       = element_text(size = 8, colour = "gray50", hjust = 0),
    axis.line          = element_line(colour = "gray30"),
    axis.ticks         = element_line(colour = "gray30"),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.y = element_line(colour = "gray90", linewidth = 0.3)
  )

print(p)

# Save
ggsave("nearest_neighbour_TauPI.pdf", plot = p, width = 7, height = 5)
ggsave("nearest_neighbour_TauPI.png", plot = p, width = 7, height = 5, dpi = 300)
message("Saved: nearest_neighbour_TauPI.pdf and .png")
