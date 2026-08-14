#!/usr/bin/env Rscript
# Tangle Burden Analysis - Publication-Ready Figures
# Violin plot + summary statistics table

library(tidyverse)
library(ggplot2)
library(ggpubr)

# ====================================================================
# LOAD DATA
# ====================================================================
# Make sure tangles.csv is in your working directory
# setwd("~/your_project_folder")

tangles <- read.csv("tangles.csv")

print("Data loaded:")
print(head(tangles))
print(paste("Total samples:", nrow(tangles)))

# ====================================================================
# SUMMARY STATISTICS TABLE
# ====================================================================
print("\n=== BURDEN STATISTICS ===")

summary_stats <- tangles %>%
  summarise(
    n = n(),
    Mean = mean(Burden_percent, na.rm = TRUE),
    SD = sd(Burden_percent, na.rm = TRUE),
    SEM = SD / sqrt(n),
    Median = median(Burden_percent, na.rm = TRUE),
    Min = min(Burden_percent, na.rm = TRUE),
    Max = max(Burden_percent, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# Save summary table
write.csv(summary_stats, "burden_summary_statistics.csv", row.names = FALSE)
print("\nSummary table saved to: burden_summary_statistics.csv")

# ====================================================================
# VIOLIN PLOT - Publication Quality
# ====================================================================

p1 <- ggplot(tangles, aes(x = "All Samples", y = Burden_percent)) +
  # Violin plot
  geom_violin(fill = "lightblue", color = "black", alpha = 0.7, linewidth = 0.8) +
  # Individual points
  geom_jitter(width = 0.1, size = 3, alpha = 0.6, color = "steelblue") +
  # Mean point
  stat_summary(fun = mean, geom = "point", size = 4, color = "red", 
               shape = "diamond", stroke = 1.5) +
  # Mean line
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, 
               color = "red", linewidth = 1) +
  # Median line
  stat_summary(fun = median, geom = "crossbar", width = 0.5, 
               color = "black", linewidth = 1, linetype = "dashed") +
  # Error bars (SEM)
  stat_summary(fun.data = function(x) {
    mean_val <- mean(x)
    se <- sd(x) / sqrt(length(x))
    data.frame(y = mean_val, ymin = mean_val - se, ymax = mean_val + se)
  }, geom = "errorbar", color = "red", width = 0.15, linewidth = 1) +
  
  # Formatting
  labs(
    title = "Tangle Burden Distribution",
    x = "",
    y = "Tangle Burden (%)",
    caption = paste("n =", nrow(tangles), "samples")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    plot.caption = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_line(color = "gray90")
  )

print(p1)

# Save plot
ggsave("tangle_burden_violin.pdf", plot = p1, width = 5, height = 7, dpi = 300)
ggsave("tangle_burden_violin.png", plot = p1, width = 5, height = 7, dpi = 300)
print("Plots saved: tangle_burden_violin.pdf and .png")

# ====================================================================
# INDIVIDUAL SAMPLE COMPARISON (if you have multiple groups)
# ====================================================================

# If your data has a "Group" column (e.g., WT, Mutant), uncomment:
# 
# p2 <- ggplot(tangles, aes(x = Group, y = Burden_percent, fill = Group)) +
#   geom_violin(alpha = 0.7) +
#   geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
#   stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
#   theme_minimal() +
#   labs(title = "Tangle Burden by Group",
#        x = "Group",
#        y = "Burden (%)") +
#   theme(legend.position = "none")
# 
# print(p2)
# ggsave("burden_by_group.pdf", plot = p2, width = 8, height = 6, dpi = 300)

# ====================================================================
# DENSITY PLOT ALTERNATIVE
# ====================================================================

p3 <- ggplot(tangles, aes(x = Burden_percent)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(Burden_percent)), 
             color = "red", linetype = "dashed", linewidth = 1,
             label = "Mean") +
  geom_vline(aes(xintercept = median(Burden_percent)), 
             color = "darkblue", linetype = "dotted", linewidth = 1,
             label = "Median") +
  labs(
    title = "Distribution of Tangle Burden",
    x = "Tangle Burden (%)",
    y = "Frequency",
    caption = paste("n =", nrow(tangles), "samples")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

print(p3)
ggsave("tangle_burden_histogram.pdf", plot = p3, width = 6, height = 5, dpi = 300)

# ====================================================================
# STATISTICAL TESTS
# ====================================================================

print("\n=== STATISTICAL TESTS ===")

# Test for normality (Shapiro-Wilk)
normality_test <- shapiro.test(tangles$Burden_percent)
print("Shapiro-Wilk normality test:")
print(normality_test)

if (normality_test$p.value > 0.05) {
  print("Data is normally distributed (p > 0.05)")
} else {
  print("Data is NOT normally distributed (p < 0.05)")
}

# ====================================================================
# CREATE PUBLICATION TABLE
# ====================================================================

# Create formatted table for paper
pub_table <- data.frame(
  Metric = c("N (samples)", "Mean ± SD (%)", "Median (%)", 
             "Range (%)", "95% CI"),
  Value = c(
    nrow(tangles),
    paste0(round(summary_stats$Mean, 4), " ± ", round(summary_stats$SD, 4)),
    round(summary_stats$Median, 4),
    paste0(round(summary_stats$Min, 4), " - ", round(summary_stats$Max, 4)),
    paste0(round(summary_stats$Mean - 1.96*summary_stats$SEM, 4), " to ",
           round(summary_stats$Mean + 1.96*summary_stats$SEM, 4))
  )
)

print("\n=== PUBLICATION TABLE ===")
print(pub_table)

write.csv(pub_table, "burden_publication_table.csv", row.names = FALSE)
print("\nPublication table saved to: burden_publication_table.csv")

print("\n=================================")
print("✓ Analysis complete!")
print("=================================")
print("Files created:")
print("  • tangle_burden_violin.pdf")
print("  • tangle_burden_violin.png")
print("  • tangle_burden_histogram.pdf")
print("  • burden_summary_statistics.csv")
print("  • burden_publication_table.csv")
print("=================================")
