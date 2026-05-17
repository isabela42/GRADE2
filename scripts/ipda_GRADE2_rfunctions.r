# ============================================================
# FUNCTION: PLOTS
# ============================================================
# Written by Isabela Almeida
# Created on May 13, 2026
# Last modified on May 17, 2026
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
# Input: plot_heat (plot_condition_col, plot_sec_condition, transcript, expr)
# Output: ggplot object
# Usage: heatmap_plot(plot_heat)
# ------------------------------------------------------------

heatmap_plot <- function(df){
  #breaks <- c(0, 0.5, 1, 3, max(df$expr, na.rm = TRUE))

  gtex_colors <- c(
    "#ffffcc",  # 0      pale yellow
    "#a0d97a",  # 0.5    light green
    "#78c679",  # 1      green
    "#41b6c4",  # 3      teal
    "#081d58"   # max    navy
  )

  #values <- scales::rescale(breaks)

  plot1 <- ggplot(plot_heat, aes(x = .data[[plot_condition]], y = .data[[plot_sec_condition]], fill = expr)) +
  geom_tile(color = "white", height = 1, width = 1, linewidth = 0.3) +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(
    colors = gtex_colors,
    #values = values,
    #limits = c(0, max(df$expr, na.rm = TRUE)),
    #oob = scales::squish,
    #breaks = breaks,
    #labels = scales::label_number(accuracy = 0.01),
    na.value = "grey92",
    guide = guide_colorbar(
      direction = "vertical",
      barwidth = 0.8,
      barheight = 10,
      ticks = TRUE,
      frame.colour = "white"
  )) +
  labs(
  fill = "log2(TPM + 1)\n(NA = grey)"
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 30, l = 5.5),
    axis.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.4, "cm"),
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
  
  plot2_df <- df %>%
  dplyr::group_by(.data[[plot_condition]]) %>%
  dplyr::summarise(
    expr = sum(expr, na.rm = TRUE),
    .groups = "drop"
  )

  plot2 <- ggplot(plot2_df, aes(x = .data[[plot_condition]], y = outstem, fill = expr)) +
  geom_tile(color = "white", height = 1, width = 1, linewidth = 0.3) +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(
    colors = gtex_colors,
    #values = values,
    #limits = c(0, max(df$expr, na.rm = TRUE)),
    #oob = scales::squish,
    #breaks = breaks,
    #labels = scales::label_number(accuracy = 0.01),
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
    legend.justification = "center",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.4, "cm"),
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

    plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 1)

  return(plot)
}