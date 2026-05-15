# ============================================================
# FUNCTION: PLOTS
# ============================================================
# Written by Isabela Almeida
# Created on May 13, 2026
# Last modified on May 14, 2026
# Version: 1.0.0
#
# DESCRIPTION: Plot functions
#
# INPUT: See each plot function for info.
#
# OUTPUT: Plot
#
# USAGE (inside script.R):
#   source("path/to/functions.R")
#   function_plot(df, script_palette)
#
# NOTES:
#   - Designed for GRADE2 pipeline use
#   - Does NOT handle file I/O or argument parsing
#   - Assumes some level of preprocessing is already done
# ============================================================

# ------------------------------------------------------------
# BOX PLOT
# Description: Log-scale distribution per group
# Input: df, plot_condition, outstem
# Output: ggplot object
# Usage: box_plot(df, plot_condition, outstem)
# ------------------------------------------------------------

box_plot <- function(df, plot_condition, outstem) {
  plot <- ggplot(df, aes(x = .data[[plot_condition]], y = .data[[outstem]])) +
  geom_boxplot(fill = "deepskyblue4", alpha=1, color = "black", size = 0.3) +
  labs(
    x = plot_condition,
    y = "log2(TPM+1)"
  ) +
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold")
  )
  return(plot)
}

# ------------------------------------------------------------
# BEESAWARM PLOT
# Description: Log-scale distribution per group
# Input: df, plot_condition, outstem
# Output: ggplot object
# Usage: beeswarm_plot(df, plot_condition, outstem)
# ------------------------------------------------------------

beeswarm_plot <- function(df, plot_condition, outstem) {
  plot <- ggplot(df, aes(x = .data[[plot_condition]], y = .data[[outstem]])) +
  geom_boxplot(fill = "deepskyblue4", alpha=1, color = "black", size = 0.3) +
  geom_quasirandom(color = "black", size = 1, alpha = 0.7, width = 0.2) +
  labs(
    x = plot_condition,
    y = "log2(TPM+1)"
  ) +
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold")
  )
  return(plot)
}

# ------------------------------------------------------------
# HEATMAP PLOT (GTEx-style)
# Description: Log-scale distribution per group
# Input: df (plot_condition_col, transcript, expr)
# Output: ggplot object
# Usage: heatmap_plot(df)
# ------------------------------------------------------------

heatmap_plot <- function(df){
  gtex_colors <- c(
    "#ffffcc",  # pale yellow
    "#c7e9b4",  # light green
    "#7fcdbb",  # teal
    "#2c7fb8",  # blue
    "#081d58"   # dark navy
  )
  
  plot <- ggplot(plot_heat, aes(x = .data[[plot_condition]], y = source, fill = expr)) +
  geom_tile(color = "white", height = 1, width = 1, linewidth = 0.3) +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(
    colors = gtex_colors,
    na.value = "grey92",
    guide = guide_colorbar(
        direction = "horizontal",
        barwidth = 10,
        barheight = 0.8,
        ticks = TRUE,
        frame.colour = "white"
  )) +
  labs(
    fill = "log2(TPM + 1)"
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 30, l = 5.5),
    axis.title = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text.y.right = element_text(
      size = 8,
      hjust = 0
      ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
    ))

  return(plot)
}