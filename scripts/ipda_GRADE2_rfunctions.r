# ============================================================
# FUNCTIONS
# ============================================================
# Written by Isabela Almeida
# Created on May 13, 2026
# Last modified on Nov 04, 2026
# Version: 1.0.0
#
# DESCRIPTION: GRADE2 R FUNCTIONS
#
# INPUT/OUTPUT: See each function for info.
#
# USAGE (inside script.R):
#   source("path/to/functions.R")
#   function(function_requirements)
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

# ------------------------------------------------------------
# BOX PLOT ALL
# Description: Log-scale distribution per group
# Input: df, plot_condition
# Output: ggplot object
# Usage: box_plot(df, plot_condition)
# ------------------------------------------------------------

box_plot_all <- function(df, plot_condition, outstem) {
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
# BOX PLOT FACET CATEGRORY
# Description: Log-scale distribution per group
# Input: df, plot_condition, outstem
# Output: ggplot object
# Usage: box_plot(df, plot_condition, outstem)
# ------------------------------------------------------------

box_plot_facet <- function(df, plot_condition, outstem) {
  plot <- ggplot(df, aes(x = .data[[plot_condition]], y = .data[[outstem]])) +
  geom_boxplot(fill = "deepskyblue4", alpha=1, color = "black", size = 0.3) +
  facet_wrap(~ category, scales = "free_x") +
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
# BEESAWARM PLOT FACET CATEGORY
# Description: Log-scale distribution per group
# Input: df, plot_condition, outstem
# Output: ggplot object
# Usage: beeswarm_plot(df, plot_condition, outstem)
# ------------------------------------------------------------

beeswarm_plot_facet <- function(df, plot_condition, outstem) {
  plot <- ggplot(df, aes(x = .data[[plot_condition]], y = .data[[outstem]])) +
  geom_boxplot(fill = "deepskyblue4", alpha=1, color = "black", size = 0.3) +
  geom_quasirandom(color = "black", size = 1, alpha = 0.7, width = 0.2) +
  facet_wrap(~ category, scales = "free_x") +
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
# HEATMAP PLOT FACET CATEGORY (GTEx-style)
# Description: Log-scale distribution per group
# Input: plot_heat (plot_condition_col, plot_sec_condition, transcript, expr)
# Output: ggplot object
# Usage: heatmap_plot(plot_heat)
# ------------------------------------------------------------

heatmap_plot_facet <- function(df){
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
  facet_wrap(~ category, scales = "free_x") +
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
  facet_wrap(~ category, scales = "free_x") +
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



################### GRADE2 071 FUNCTIONS
annotate_qlfo <- function(qlft, anno, level) {
  qlfo <- cbind(feature = rownames(qlft), qlft)
  rownames(qlfo) <- rownames(qlft)
  if (level == "transcript") {
    qlfo <- left_join(qlfo, anno[, c("TRANSCRIPT", "GENE", "SYMBOLGENE", "SYMBOLTRANSCRIPT")], by = c("feature" = "TRANSCRIPT"))
    qlfo$SYMBOL <- qlfo$SYMBOLGENE
  } else {qlfo <- left_join(qlfo, anno[, c("GENE", "SYMBOLGENE")], by = c("feature" = "GENE"))}
  qlfo$diffexpressed <- "NO"
  qlfo$diffexpressed[qlfo$logFC >  lfc_cut & qlfo$FDR < fdr_cut] <- "UP"
  qlfo$diffexpressed[qlfo$logFC < -lfc_cut & qlfo$FDR < fdr_cut] <- "DOWN"
  qlfo$pathway <- NA_character_
  for (pw in names(dna_repair_genes)) {qlfo$pathway[!is.na(qlfo$SYMBOL) & qlfo$SYMBOL %in% dna_repair_genes[[pw]]] <- pw}
  qlfo
}

annotate_qlfo_ensemblr <- function(qlft, anno) {
  qlfo           <- cbind(feature=rownames(qlft), qlft)
  rownames(qlfo) <- rownames(qlft)
  qlfo$FeatClean <- sub("\\.[0-9]+$", "", qlfo$feature)
  qlfo <- left_join(qlfo, anno[, c("ENSEMBL","SYMBOL","ENTREZID")],
                    by=c("FeatClean"="ENSEMBL"))
  qlfo <- qlfo %>% relocate(SYMBOL, ENTREZID, .after=feature)
  qlfo$diffexpressed <- "NO"
  qlfo$diffexpressed[qlfo$logFC >  lfc_cut & qlfo$FDR < fdr_cut] <- "UP"
  qlfo$diffexpressed[qlfo$logFC < -lfc_cut & qlfo$FDR < fdr_cut] <- "DOWN"
  # Pathway annotation
  qlfo$pathway <- NA_character_
  for (pw in names(dna_repair_genes)) {
    qlfo$pathway[!is.na(qlfo$SYMBOL) & qlfo$SYMBOL %in% dna_repair_genes[[pw]]] <- pw
  }
  qlfo
}

plot_pca <- function(y, grps, ctr, outstem, outdir, level) {
  log2cpm <- as.data.frame(edgeR::cpm(y, normalized.lib.sizes=TRUE,
                                       log=TRUE, prior.count=2))
  rv  <- matrixStats::rowVars(as.matrix(log2cpm))
  sel <- rownames(log2cpm)[order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]]

  pca_res <- prcomp(t(log2cpm[sel, ]), center=TRUE, scale.=FALSE)
  pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  pca_df  <- as.data.frame(pca_res$x[, 1:min(2, ncol(pca_res$x))])
  pca_df$sample <- rownames(pca_df)
  pca_df$group  <- grps

  lvls    <- levels(factor(grps))
  pal     <- RColorBrewer::brewer.pal(max(3, length(lvls)), "Set1")[seq_along(lvls)]
  col_key <- setNames(pal, lvls)

  p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group, label=sample)) +
    # Dotted crosshair lines per point
    geom_segment(aes(x=PC1, xend=PC1, y=min(.data$PC2)-5, yend=PC2),
                 linetype="dotted", alpha=0.3, color="grey60", linewidth=0.3) +
    geom_segment(aes(x=min(.data$PC1)-5, xend=PC1, y=PC2, yend=PC2),
                 linetype="dotted", alpha=0.3, color="grey60", linewidth=0.3) +
    geom_point(size=4, alpha=0.9) +
    geom_text_repel(size=3, show.legend=FALSE,
                    box.padding=0.4, point.padding=0.3) +
    scale_color_manual(values=col_key) +
    labs(x=paste0("PC1 (", pct_var[1], "%)"),
         y=paste0("PC2 (", pct_var[2], "%)"),
         title=paste("PCA -", outstem), color="Group") +
    theme_classic() +
    theme(legend.position="right",
          panel.grid.major=element_line(color="grey92", linewidth=0.3))

  scree_df <- data.frame(
    PC=paste0("PC", seq_along(pct_var)),
    Variance=pct_var
  )[seq_len(min(10, length(pct_var))), ]
  scree_df$PC <- factor(scree_df$PC, levels=scree_df$PC)

  p_scree <- ggplot(scree_df, aes(x=PC, y=Variance)) +
    geom_col(fill="steelblue") +
    labs(x=NULL, y="% Variance explained", title="Scree plot") +
    theme_classic()

  pdf(file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".pca.pdf")), width=7, height=6)
    print(p_pca)
    print(p_scree)
  dev.off()
  cat("## PCA saved:", stem, "\n")
}


plot_volcano <- function(qlfo, ctr, outstem, outdir, label_mode, level) {
  n_up   <- sum(qlfo$diffexpressed=="UP",   na.rm=TRUE)
  n_down <- sum(qlfo$diffexpressed=="DOWN",  na.rm=TRUE)

  # Top DEG labels (top 5, excluding pathway genes to avoid overlap)
  top_n <- 5
  label_df <- bind_rows(
    qlfo %>% filter(diffexpressed=="UP",   is.na(pathway)) %>%
      arrange(FDR) %>% slice_head(n=top_n),
    qlfo %>% filter(diffexpressed=="DOWN", is.na(pathway)) %>%
      arrange(FDR) %>% slice_head(n=top_n)
  )

  # Pathway DEG genes
  path_deg <- qlfo %>% filter(!is.na(pathway) & diffexpressed != "NO")

  base_colors <- c("DOWN"="#00AFBB", "NO"="grey80", "UP"="#bb0c00")

  p <- ggplot(qlfo, aes(x=logFC, y=-log10(FDR))) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), col="gray", linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut),       col="gray", linetype="dashed") +
    # All points colored by DE status
    geom_point(aes(color=diffexpressed), size=1.5, alpha=0.5) +
    scale_color_manual(name="DEG", values=base_colors,
      labels=c("Downregulated","Not significant","Upregulated")) +
    # Pathway DEG points — outlined on top, colored by pathway
    ggnewscale::new_scale_color() +
    geom_point(data=path_deg, aes(color=pathway),
               size=3, shape=21, stroke=1.2, fill=NA) +
    scale_color_manual(name="Pathway", values=pathway_colors) +
    # Top DEG labels — black
    geom_text_repel(data=label_df, aes(label=SYMBOL), color="black",
      size=3, fontface="plain", box.padding=0.5, point.padding=0.5,
      segment.color="grey50", max.overlaps=20, show.legend=FALSE) +
    # Pathway DEG labels — colored bold
    geom_text_repel(data=path_deg, aes(label=SYMBOL, color=pathway),
      size=3, fontface="bold", box.padding=0.5, point.padding=0.5,
      segment.color="grey50", segment.linetype="dashed",
      max.overlaps=20, show.legend=FALSE) +
    coord_cartesian(xlim=c(-10, 10)) +
    scale_x_continuous(breaks=seq(-10, 10, 2)) +
    labs(x=expression("log"[2]*"FC"), y=expression("-log"[10]*"FDR"),
         title=paste("Volcano -", outstem),
         subtitle=paste0("UP: ", n_up, "  |  DOWN: ", n_down,
                         "  |  FDR<", fdr_cut, " & |logFC|>", lfc_cut)) +
    theme_classic() +
    theme(legend.position="right")

  pdf(file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".volcano.pdf")), width=9, height=6)
    print(p)
  dev.off()
  cat("## Volcano saved:", stem, "\n")
}

run_fgsea <- function(qlfo, ctr, outstem, input_gmt, outdir, level) {
  log_vec        <- qlfo$logFC
  names(log_vec) <- qlfo$SYMBOL
  log_vec        <- log_vec[!is.na(names(log_vec)) & is.finite(log_vec)]
  log_vec        <- log_vec[!duplicated(names(log_vec))]

  gmt_files <- list(
    KEGG     = file.path(input_gmt, "c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt"),
    Reactome = file.path(input_gmt, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"),
    GO_BP    = file.path(input_gmt, "c5.go.bp.v2025.1.Hs.symbols.gmt")
  )

  for (db in names(gmt_files)) {
    if (!file.exists(gmt_files[[db]])) {
      cat("## WARNING: GMT not found:", gmt_files[[db]], "\n"); next
    }
    cat("## fgsea:", db, "...\n")
    pathways <- fgsea::gmtPathways(gmt_files[[db]])
    set.seed(42)
    res <- fgsea::fgsea(pathways=pathways, stats=log_vec, minSize=15, maxSize=500)
    res <- res[order(res$padj), ]
    write.table(as.data.frame(res) %>% dplyr::select(-leadingEdge), file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".", db, ".fgsea.tsv")),
      quote=FALSE, row.names=FALSE, sep="\t")

    res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
    if (nrow(res_sig) == 0) { cat("## No sig pathways:", db, "\n"); next }

    res_sig_df <- as.data.frame(res_sig)[, !colnames(res_sig) %in% "leadingEdge"]
    plot_df <- unique(rbind(
      head(res_sig_df[order(-res_sig_df$NES), ], 10),
      head(res_sig_df[order( res_sig_df$NES), ], 10)
    ))
    plot_df$pathway <- factor(plot_df$pathway,
                               levels=plot_df$pathway[order(plot_df$NES)])
    # Clean pathway names for display
    plot_df$label <- gsub(paste0("^", db, "_"), "", plot_df$pathway)
    plot_df$label <- gsub("_", " ", plot_df$label)
    plot_df$label <- factor(plot_df$label, levels=plot_df$label[order(plot_df$NES)])

    # Dotplot: size = number of genes in pathway, color = NES direction
    p_path <- ggplot(plot_df, aes(x=NES, y=label)) +
      geom_point(aes(size=size, color=NES), alpha=0.85) +
      scale_color_gradient2(
        low="#2166AC", mid="white", high="#BB0C00",
        midpoint=0, name="NES"
      ) +
      scale_size_continuous(name="Gene set size", range=c(3, 10)) +
      geom_vline(xintercept=0, linetype="dashed", color="grey50") +
      labs(x="Normalized Enrichment Score", y=NULL,
           title=paste(db, "-", outstem),
           subtitle=paste0(outstem, "\npadj<0.05, top 10 up + top 10 down")) +
      theme_minimal(base_size=11) +
      theme(axis.text.y=element_text(size=8),
            panel.grid.major.y=element_line(color="grey92"),
            panel.grid.major.x=element_line(color="grey85"))

    h <- max(6, nrow(plot_df) * 0.38 + 2)
    pdf(file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".", db, ".fgsea.pdf")), width=11, height=h)
      print(p_path)
    dev.off()
    cat("##", db, "dotplot saved\n")
  }
}

plot_heatmap_de <- function(tpmcnt, groups_vec, genesets, qlfo_list,
                         outstem, level, outdir,annotpm, substem) {
  if (is.null(genesets)) return(invisible(NULL))

  # Prepare TPM matrix
  tpm_df <- tpmcnt
  tpm_df$feature <- rownames(tpm_df)
  if (level == "transcript") {
    tpm_df <- left_join(tpm_df, annotpm[, c("TRANSCRIPT", "GENE", "SYMBOLGENE", "SYMBOLTRANSCRIPT")], by = c("feature" = "TRANSCRIPT"))
    tpm_df$SYMBOL <- tpm_df$SYMBOLGENE
    tpm_df <- tpm_df %>% filter(!is.na(SYMBOL)) %>% distinct(SYMBOL, .keep_all=TRUE)
    rownames(tpm_df) <- tpm_df$SYMBOL
    expr_cols <- colnames(tpmcnt)    
  } else {
    tpm_df <- left_join(tpm_df, annotpm[, c("GENE", "SYMBOLGENE")], by = c("feature" = "GENE"))
    tpm_df$SYMBOL <- tpm_df$SYMBOLGENE
    tpm_df <- tpm_df %>% filter(!is.na(SYMBOL)) %>% distinct(SYMBOL, .keep_all=TRUE)
    rownames(tpm_df) <- tpm_df$SYMBOL
    expr_cols <- colnames(tpmcnt)
  }
  tpm_mat <- as.matrix(tpm_df[, expr_cols])
  storage.mode(tpm_mat) <- "numeric"

  use_top10 <- "top10DEG" %in% genesets
  path_sets <- genesets[genesets != "top10DEG"]
  if ("ALL" %in% path_sets) path_sets <- names(dna_repair_genes)

  # Column annotation — ASO1 and ASO2 same color
  cond_lvls <- unique(groups_vec)
  # Auto-detect KO groups (non-CTRL) and merge their color
  ctrl_grp  <- grep("CTRL|ctrl|NC|nc|NEG|neg", cond_lvls, value=TRUE)
  ko_grps   <- setdiff(cond_lvls, ctrl_grp)
  # Assign colors: CTRL=grey, KO groups all same red
  cond_cols <- setNames(
    c(rep("#BB0C00", length(ko_grps)), rep("#666666", length(ctrl_grp))),
    c(ko_grps, ctrl_grp)
  )
  col_ann <- data.frame(condition=groups_vec, row.names=colnames(tpm_mat))

  make_heatmap <- function(mat, row_ann, stem_suffix, path_lvls, substem, level) {
    zv  <- apply(mat, 1, var) == 0
    mat <- mat[!zv, , drop=FALSE]
    if (nrow(mat) < 2) { cat("## WARNING: Too few rows for heatmap", stem_suffix, "\n"); return() }

    # Row annotation aligned
    if (!is.null(row_ann)) {
      row_ann <- row_ann[rownames(mat), , drop=FALSE]
    }

    path_pal  <- c(pathway_colors, top10DEG="#aaaaaa")
    used_lvls <- if (!is.null(row_ann)) unique(row_ann$pathway) else character(0)
    p_cols    <- path_pal[used_lvls]
    p_cols    <- p_cols[!is.na(p_cols)]

    ann_colors <- list(condition=cond_cols)
    if (length(p_cols) > 0) ann_colors$pathway <- p_cols

    h    <- max(8, nrow(mat) * 0.28 + 3)
    pdf(file.path(outdir, paste0(outstem, ".", level, ".", substem, ".", stem_suffix, ".heatmap.pdf")), width=11, height=h)
      pheatmap(
        mat                      = mat,
        scale                    = "row",
        show_rownames            = TRUE,
        show_colnames            = TRUE,
        cluster_rows             = TRUE,
        cluster_cols             = TRUE,
        clustering_method        = "complete",
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "euclidean",
        annotation_row           = row_ann,
        annotation_col           = col_ann,
        annotation_colors        = ann_colors,
        cutree_rows              = max(1, length(path_lvls)),
        cutree_cols              = length(cond_lvls),
        fontsize_row             = 10,
        fontsize_col             = 9,
        main                     = paste("Heatmap -", outstem)
      )
    dev.off()
  }

  # ── Heatmap 1: pathway genes ──────────────────────────────────────────────
  if (length(path_sets) > 0) {
    path_df <- bind_rows(lapply(path_sets, function(p) {
      if (!p %in% names(dna_repair_genes)) return(NULL)
      data.frame(SYMBOL=dna_repair_genes[[p]], pathway=p, stringsAsFactors=FALSE)
    })) %>% distinct(SYMBOL, .keep_all=TRUE)

    genes_found <- intersect(path_df$SYMBOL, rownames(tpm_mat))
    cat("## Pathway heatmap genes found:", length(genes_found), "/", nrow(path_df), "\n")

    if (length(genes_found) >= 2) {
      mat     <- tpm_mat[genes_found, , drop=FALSE]
      row_ann <- path_df %>%
        filter(SYMBOL %in% rownames(mat)) %>%
        distinct(SYMBOL, .keep_all=TRUE) %>%
        column_to_rownames("SYMBOL") %>%
        dplyr::select(pathway)
      make_heatmap(mat, row_ann,
                   paste0(paste(path_sets, collapse="-"), "_pathways"),
                   path_sets, substem, level)
    }
  }

  # ── Heatmap 2: top 10 DEG ─────────────────────────────────────────────────
  if (use_top10 && length(qlfo_list) > 0) {
    top_deg <- bind_rows(lapply(qlfo_list, function(q) {
      bind_rows(
        q %>% filter(diffexpressed=="UP")   %>% arrange(FDR) %>% slice_head(n=10),
        q %>% filter(diffexpressed=="DOWN") %>% arrange(FDR) %>% slice_head(n=10)
      )
    })) %>% filter(!is.na(SYMBOL)) %>% distinct(SYMBOL) %>% pull(SYMBOL)

    genes_found <- intersect(top_deg, rownames(tpm_mat))
    cat("## Top10DEG heatmap genes found:", length(genes_found), "/", length(top_deg), "\n")

    if (length(genes_found) >= 2) {
      genes_found  <- unique(genes_found)  # ensure no duplicates
      mat          <- tpm_mat[genes_found[genes_found %in% rownames(tpm_mat)], , drop=FALSE]
      row_ann      <- data.frame(pathway=rep("top10DEG", nrow(mat)),
                                 row.names=rownames(mat))
      make_heatmap(mat, row_ann, "top10DEG", "top10DEG", substem, level)
    }
  }
}

run_edger <- function(counts_mat, groups, contrast_str) {
  parts  <- strsplit(contrast_str, "_vs_")[[1]]
  grp_ko <- parts[1]; grp_ref <- parts[2]
  keep_s <- groups %in% c(grp_ko, grp_ref)
  mat    <- counts_mat[, keep_s, drop=FALSE]
  grps   <- factor(groups[keep_s], levels=c(grp_ref, grp_ko))

  cat("\n## Contrast:", contrast_str,
      "| n(", grp_ko, ")=", sum(grps==grp_ko),
      "vs n(", grp_ref, ")=", sum(grps==grp_ref), "\n")

  y      <- DGEList(counts=mat, group=grps)
  keep   <- rowSums(cpm(y) > min_cpm) >= min_samples
  y      <- y[keep, , keep.lib.sizes=FALSE]
  cat("## Features after filter:", nrow(y), "\n")

  y      <- calcNormFactors(y)
  design <- model.matrix(~grps)
  y      <- estimateGLMCommonDisp(y, design)
  y      <- estimateGLMTrendedDisp(y, design)
  y      <- estimateGLMTagwiseDisp(y, design)
  fit    <- glmQLFit(y, design)
  qlf    <- glmQLFTest(fit, coef=2)

  nr     <- nrow(qlf)
  qlft   <- topTags(qlf, adjust.method="BH", n=nr)$table
  p      <- qlft[, ncol(qlft)-1]
  qlft$AdjPValue <- p.adjust(p, method="BH", n=ngenes)

  list(qlft=qlft, y=y)
}

## 083 plots
## ---------------------------------------------------------------- utilities
# Split a comma-separated string into a clean character vector ("" -> character(0))
split_terms <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(character(0))
  x <- trimws(strsplit(x, ",")[[1]])
  x[nzchar(x)]
}
 
# Open a PDF, run code, and ALWAYS close the device (even on error)
with_pdf <- function(file, width, height, code) {
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(code)
  invisible(file)
}
 
# Drop metadata rows whose any column matches one of the terms
exclude_samples <- function(metadata, terms) {
  if (length(terms) == 0) return(metadata)
  row_txt <- do.call(paste, c(metadata, sep = " | "))
  drop <- grepl(paste(terms, collapse = "|"), row_txt, ignore.case = TRUE)
  message("Excluding ", sum(drop), " samples matching: ", paste(terms, collapse = ", "))
  metadata[!drop, , drop = FALSE]
}
 
## ---------------------------------------------------------------- STAR logs
read_star_log <- function(log_file){
  cols <- c("input_reads", "read_length", "unique_reads", "unique_map_pct", "multi_reads", "multimapped_pct",
            "too_many_loci_reads", "too_many_loci_pct", "unmapped_mismatch_reads", "unmapped_mismatch_pct",
            "unmapped_short_reads", "unmapped_short_pct", "unmapped_other_reads", "unmapped_other_pct",
            "chimeric_reads", "chimeric_pct")
  if (!file.exists(log_file)) {
    warning("STAR log not found: ", log_file)
    return(as.data.frame(setNames(as.list(rep(NA_real_, length(cols))), cols)))
  }
  x <- readLines(log_file)
  get_val <- function(pattern){
    hit <- grep(pattern, x, value = TRUE)
    if(length(hit) == 0) {return(NA)}
    trimws(sub(".*\\|", "", hit[1]))
  }
  pct <- function(pattern){as.numeric(gsub("%", "", get_val(pattern)))}
  num <- function(pattern){as.numeric(gsub(",", "", get_val(pattern)))}
  data.frame(
    input_reads              = num("Number of input reads"),
    read_length              = num("Average input read length"),
    unique_reads             = num("Uniquely mapped reads number"),
    unique_map_pct           = pct("Uniquely mapped reads %"),
    multi_reads              = num("Number of reads mapped to multiple loci"),
    multimapped_pct          = pct("% of reads mapped to multiple loci"),
    too_many_loci_reads      = num("Number of reads mapped to too many loci"),
    too_many_loci_pct        = pct("% of reads mapped to too many loci"),
    unmapped_mismatch_reads  = num("Number of reads unmapped: too many mismatches"),
    unmapped_mismatch_pct    = pct("% of reads unmapped: too many mismatches"),
    unmapped_short_reads     = num("Number of reads unmapped: too short"),
    unmapped_short_pct       = pct("% of reads unmapped: too short"),
    unmapped_other_reads     = num("Number of reads unmapped: other"),
    unmapped_other_pct       = pct("% of reads unmapped: other"),
    chimeric_reads           = num("Number of chimeric reads"),
    chimeric_pct             = pct("% of chimeric reads")
  )
}
 
# Table of STAR metrics per sample (one row per metadata row)
build_star_qc <- function(metadata, star_dir) {
  metadata %>%
    mutate(log_file = file.path(star_dir, paste0("salign-star-index_", sequencing.id),
                                paste0("salign-star-index_", sequencing.id, "Log.final.out"))) %>%
    rowwise() %>% mutate(qc = list(read_star_log(log_file))) %>% ungroup() %>% tidyr::unnest(qc)
}
 
## ---------------------------------------------------------------- themes / plots
my_theme <- function(axis_x_size = 8, axis_y_size = 8, axis_title_size = 12, strip_size = 12, legend_size = 12, panel_spacing = 1){
    theme_bw() +
    theme(
        panel.spacing = unit(panel_spacing, "lines"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_x_size),
        axis.text.y = element_text(size = axis_y_size),
        axis.title.x = element_text(size = axis_title_size),
        axis.title.y = element_text(size = axis_title_size),
        plot.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size),
        strip.text = element_text(size = strip_size, face = "bold")
    )
}
 
metadata_barplot_facet <- function(df, x_condition, y_condition, fill_condition = NULL) {
  if (is.null(fill_condition)) {
    plot_df <- df %>% count(.data[[y_condition]], .data[[x_condition]]) %>%
      mutate(x_ordered = tidytext::reorder_within(.data[[x_condition]], -n, .data[[y_condition]]))
    plot <- ggplot(plot_df, aes(x = x_ordered, y = n)) +
    geom_col(fill = "#0072B2", alpha=1, color = "black") +
    geom_text(aes(label = n), vjust = -0.2, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    tidytext::scale_x_reordered(expand = c(0, 0))
  } else {
    plot_df <- df %>% count(.data[[y_condition]], .data[[x_condition]], .data[[fill_condition]]) %>%
      group_by(.data[[y_condition]], .data[[x_condition]]) %>%
      mutate(total_n = sum(n), n_segments = n()) %>%
      ungroup() %>% mutate(x_ordered = factor(.data[[x_condition]], levels = sort(unique(df[[x_condition]]))))
    plot <- ggplot(plot_df, aes(x = x_ordered, y = n, fill = .data[[fill_condition]])) +
    geom_col(position = "fill", color = "black", linewidth = 0.2) +
    geom_text(aes(label = ifelse(n > 10 | n_segments == 1, n, "")), position = position_fill(vjust = 0.5), size = 2.5, color = "white") +
    scale_y_continuous(labels = scales::percent)
  }
  plot <- plot +
    facet_grid(cols = vars(.data[[y_condition]]), scales = "free_x", space = "free_x") +
    paletteer::scale_fill_paletteer_d("ggsci::category20b_d3") +
    labs(x = x_condition, y = ifelse(is.null(fill_condition), "Samples", "Sample proportion"), fill = fill_condition) +
    theme_bw() +
    theme(
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8, margin = margin(r = 0)),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold")
    )
  return(plot)
}
 
histogram_plot <- function(df, condition, wrap = NULL, title, scales = NULL, log10 = FALSE, bins = 50){
    plot <- ggplot(df, aes(x = .data[[condition]])) +
    geom_histogram(fill = "#0072B2", bins = bins)
    if(log10){plot <- plot + scale_x_log10()}
    if(!is.null(wrap)){plot <- plot + facet_wrap(wrap, scales = scales)}
    plot + ggtitle(title) + my_theme()
}
 
barplot_distribution <- function(df, conditionX, conditionY = NULL, wrap = NULL, title, scales = "fixed"){
    if(is.null(conditionY)){
        plot <- ggplot(df, aes(x = .data[[conditionX]])) + geom_bar(fill = "#0072B2")
    } else {
        plot <- ggplot(df, aes(x = .data[[conditionX]], y = .data[[conditionY]])) + geom_col(fill = "#0072B2")}
    if(!is.null(wrap)){plot <- plot + facet_wrap(wrap, scales = scales)}
    plot + ggtitle(title) + my_theme()
}
 
boxplot_distribution <- function(df, conditionX, conditionY = NULL, wrap = NULL, title, scales = "fixed"){
    if(is.null(conditionY)){
        plot <- ggplot(df, aes(x = .data[[conditionX]])) + geom_boxplot(fill = "#0072B2")
    } else {
        plot <- ggplot(df, aes(x = .data[[conditionX]], y = .data[[conditionY]])) + geom_boxplot(fill = "#0072B2")}
    if(!is.null(wrap)){plot <- plot + facet_wrap(wrap, scales = scales)}
    plot + ggtitle(title) + my_theme()
}
 
scatter_plot <- function(df, conditionX, conditionY, gridcol = NULL, gridrow = NULL, shape = NULL, title, log10_x = FALSE){
    plot <- ggplot(df, aes(x = .data[[conditionX]], y = .data[[conditionY]]))
    if (!is.null(shape)) {plot <- plot +
        geom_point(aes(shape = .data[[shape]]), color = "#0072B2", alpha = 0.5, size = 2) +
        scale_shape_manual(values = rep_len(0:14, dplyr::n_distinct(df[[shape]])))   # recycles if > 15 levels
    } else {plot <- plot + geom_point(color = "#0072B2", alpha = 0.5)}
    if(log10_x){plot <- plot + scale_x_log10()}
    if(!is.null(gridrow) && !is.null(gridcol)){plot <- plot +
        facet_grid(rows = vars(.data[[gridrow]]), cols = vars(.data[[gridcol]]))
    } else if(!is.null(gridcol)){plot <- plot + facet_grid(cols = vars(.data[[gridcol]]))}
    plot + ggtitle(title) + my_theme()
}
 
violin_distribution <- function(df, conditionX, conditionY, gridrow = NULL, gridcol = NULL, title, scales = "fixed", show_n = FALSE,
                                plot_width = 11.67, max_n_size = 3, min_n_size = 1.2, max_strip_size = 12, min_strip_size = 4) {
  plot <- ggplot(df, aes(x = .data[[conditionX]], y = .data[[conditionY]])) +
    geom_violin(fill = "#0072B2", alpha = 0.7, color = "black", scale = "width") +
    geom_boxplot(width = 0.1) +
    ggbeeswarm::geom_quasirandom(color = "black", size = 0.5, alpha = 0.2, width = 0.2)
 
  # --- Estimate layout: violins per column panel, and inches available per violin / per panel
  free_x <- !is.null(gridcol) && scales %in% c("free", "free_x")
  n_x_all <- length(unique(df[[conditionX]]))
  if (is.null(gridcol)) {
    xk <- n_x_all
    col_names <- NULL
  } else {
    col_vals <- droplevels(factor(df[[gridcol]]))
    col_names <- levels(col_vals)
    xk <- vapply(col_names, function(cn) length(unique(df[[conditionX]][which(col_vals == cn)])), integer(1))
    if (!free_x) xk[] <- n_x_all
  }
  usable_w <- max(plot_width - 1.2 - ifelse(is.null(gridrow), 0, 0.5) - 0.2 * (length(xk) - 1), 1)
  slot_w  <- usable_w / sum(xk + 1.2)   # inches per violin
  panel_w <- slot_w * (xk + 1.2)        # inches per column panel
 
  # --- Column strip titles sized to fit their panel
  use_h4x <- !is.null(gridcol) && requireNamespace("ggh4x", quietly = TRUE)
  if (!is.null(gridcol)) {
    strip_sizes <- pmin(max_strip_size, pmax(min_strip_size, 0.9 * panel_w * 72.27 / (nchar(col_names) * 0.62)))
  }
  facet_fun <- if (use_h4x) {
    function(...) ggh4x::facet_grid2(..., strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(size = strip_sizes)))
  } else {
    facet_grid
  }
 
  if (!is.null(gridrow) && !is.null(gridcol)) {
    plot <- plot + facet_fun(rows = vars(.data[[gridrow]]), cols = vars(.data[[gridcol]]), scales = scales, space = scales)
  } else if (!is.null(gridrow)) {
    plot <- plot + facet_grid(rows = vars(.data[[gridrow]]))
  } else if (!is.null(gridcol)) {
    plot <- plot + facet_fun(cols = vars(.data[[gridcol]]), scales = scales, space = scales)
  }
 
  # --- n labels: skipped when mean is 0, size scaled to the space available per violin
  if (show_n) {
    group_vars <- c(conditionX, gridrow, gridcol)
    counts_df <- df %>% group_by(across(all_of(group_vars))) %>%
      summarise(n = n(),
                mean_y = mean(.data[[conditionY]], na.rm = TRUE),
                ymax = suppressWarnings(max(.data[[conditionY]], na.rm = TRUE)),
                .groups = "drop") %>%
      filter(is.finite(mean_y)) # To hide when mean is 0: filter(is.finite(mean_y), mean_y != 0)
    if (nrow(counts_df) > 0) {
      n_size <- min(max_n_size, max(min_n_size, 40 * slot_w / max(nchar(counts_df$n))))
      plot <- plot +
        geom_text(data = counts_df, aes(x = .data[[conditionX]], y = ymax, label = n),
                  inherit.aes = FALSE, size = n_size, vjust = -1, check_overlap = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    }
  }
 
  plot <- plot + labs(x = conditionX, y = title) + my_theme()
  # Fallback if ggh4x is not installed: one strip size for all columns (the smallest needed)
  if (!is.null(gridcol) && !use_h4x) plot <- plot + theme(strip.text.x = element_text(size = min(strip_sizes)))
  return(plot)
}
 
## ---------------------------------------------------------------- plot_df preparation
# plot_group = phenotype, replaced by phenotype.subtype for rows matching a highlight term;
# highlighted groups come first in the factor levels (in term order).
add_plot_group <- function(df, highlight_terms, type_levels) {
  df$plot_group <- df$phenotype
  if (length(highlight_terms) > 0) {
    hits <- grepl(paste(highlight_terms, collapse = "|"), df$phenotype, ignore.case = TRUE)
    df$plot_group[hits] <- df$phenotype.subtype[hits]
  }
  highlight_levels <- c()
  for (term in highlight_terms) {
    tmp <- df %>% filter(grepl(term, phenotype, ignore.case = TRUE)) %>% pull(plot_group) %>% unique() %>% sort()
    highlight_levels <- c(highlight_levels, tmp)
  }
  other_levels <- setdiff(sort(unique(df$plot_group)), highlight_levels)
  df$plot_group <- factor(df$plot_group, levels = unique(c(highlight_levels, other_levels)))
  # Unknown types are appended instead of silently becoming NA
  df$type <- factor(as.character(df$type), levels = union(type_levels, unique(as.character(df$type))))
  df$type <- droplevels(df$type)
  df
}
 
prep_identity <- function(df, terms, treat_col = "treatment") {
  df <- df[df$type == "cell-line", ]
  tr <- as.character(df[[treat_col]])
  df$treatment_status <- factor(ifelse(is.na(tr) | trimws(tr) %in% c("", "NA"), "untreated", "treated"),
                                levels = c("untreated", "treated"))
  df$identity <- as.character(df$phenotype.identity)
  hl <- c()
  for (term in terms) hl <- c(hl, sort(unique(df$identity[grepl(term, df$phenotype, ignore.case = TRUE)])))
  hl <- unique(hl)
  df$identity <- factor(df$identity, levels = c(hl, setdiff(sort(unique(df$identity)), hl)))
  df
}
 
## ---------------------------------------------------------------- STAR QC pdf
# Each plot gets its own page height; pages are joined with pdfunite (poppler-utils)
save_plot_list_pdf <- function(plots, heights, file, width = 8.27) {
  tmpdir <- tempfile("star_pages_"); dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  pdfs <- file.path(tmpdir, sprintf("%02d_%s.pdf", seq_along(plots), names(plots)))
  for (i in seq_along(plots)) ggsave(filename = pdfs[i], plot = plots[[i]], width = width, height = heights[[names(plots)[i]]])
  if (nzchar(Sys.which("pdfunite"))) {
    status <- system2("pdfunite", c(shQuote(pdfs), shQuote(file)))
    if (status != 0) warning("pdfunite failed for ", file)
  } else {
    warning("pdfunite not found in PATH; writing individual pages next to ", file)
    file.copy(pdfs, dirname(file), overwrite = TRUE)
  }
  invisible(file)
}
 
## ---------------------------------------------------------------- expression plots
plot_feature_pdf <- function(feature, file_stem, plot_df, outdir, outstem, types, highlight_terms) {
  df <- plot_df
  df$feature_value <- df[[feature]]
  ylab <- "log2(tpm + 1)"
  with_pdf(paste0(outdir, "/", outstem, ".violin.", file_stem, ".pdf"), width = 11.67, height = 5, {
    print(violin_distribution(df, conditionX = "plot_group", conditionY = "feature_value", gridcol = "category", scales = "free_x", title = ylab, show_n = TRUE) +
      ggtitle(paste(feature, "- all types combined")))
    id_df <- prep_identity(df, highlight_terms)
    if (nrow(id_df) > 0) print(violin_distribution(id_df, conditionX = "identity", conditionY = "feature_value", gridrow = "treatment_status", title = ylab, show_n = TRUE, plot_width = 24) +
      ggtitle(paste(feature, "- cell-line names")))
    for (t in types) {
      sub_df <- df[df$type == t, ]
      if (nrow(sub_df) == 0) next
      print(violin_distribution(sub_df, conditionX = "plot_group", conditionY = "feature_value", gridcol = "category", gridrow = NULL, scales = "free_x", title = ylab, show_n = TRUE) +
        ggtitle(paste(feature, "-", as.character(t))))
    }
  })
}
 
plot_mean_pdf <- function(plot_df, feature_labels, file, types, highlight_terms) {
  plot_df$feature_mean <- rowMeans(plot_df[, feature_labels, drop = FALSE], na.rm = TRUE)
  ttl <- "mean log2(tpm + 1) expression"
  with_pdf(file, width = 11.67, height = 5, {
    print(violin_distribution(plot_df, conditionX = "plot_group", conditionY = "feature_mean", gridcol = "category", scales = "free_x", title = ttl, show_n = TRUE) +
      ggtitle("All types combined"))
    id_df <- prep_identity(plot_df, highlight_terms)
    if (nrow(id_df) > 0) print(violin_distribution(id_df, conditionX = "identity", conditionY = "feature_mean", gridrow = "treatment_status", title = ttl, show_n = TRUE, plot_width = 24) +
      ggtitle("cell-line names"))
    for (t in types) {
      sub_df <- plot_df[plot_df$type == t, ]
      if (nrow(sub_df) == 0) next
      print(violin_distribution(sub_df, conditionX = "plot_group", conditionY = "feature_mean", gridcol = "category", gridrow = NULL, scales = "free_x", title = ttl, show_n = TRUE) +
        ggtitle(as.character(t)))
    }
  })
}
 
## ---------------------------------------------------------------- heatmap
heatmap_palettes <- list(
  # Sequential (use with scale = "none", i.e. absolute log2(tpm + 1))
  blues  = colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(100),
  gtex   = colorRampPalette(c("#FFFFFF", "#DCE9F5", "#7FB2DA", "#2E75B6", "#0B2C5A"))(100),
  mako   = viridis::mako(100, direction = -1),
  ylgnbu = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
  # Diverging (use with scale = "row", i.e. z-scores centred on 0)
  div      = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  div_blue = colorRampPalette(c("#08306B", "#4292C6", "#F7FBFF", "#FDD49E", "#B30000"))(100),
  div_strong = colorRampPalette(c("#053061", "#2166AC", "#67A9CF", "#F7F7F7", "#FDB863", "#E66101", "#7F3B08"))(100)
)
 
plot_heatmap <- function(plot_df, feature_labels, file, pal = heatmap_palettes$div_strong) {
  grp <- plot_df %>%
    group_by(category, plot_group) %>%
    summarise(n_samples = n(), across(all_of(feature_labels), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  heat_mat <- t(as.matrix(grp[, feature_labels]))           # rows = features, cols = groups
  colnames(heat_mat) <- make.unique(as.character(grp$plot_group))
  anno <- data.frame(category = grp$category, n_samples = grp$n_samples, row.names = colnames(heat_mat))
  heat_mat <- heat_mat[apply(heat_mat, 1, var, na.rm = TRUE) > 0, , drop = FALSE]  # no-variance features can't be row-scaled/clustered
  if (nrow(heat_mat) < 2 || ncol(heat_mat) < 2) {
    message("Heatmap skipped: need >= 2 variable features and >= 2 groups (have ", nrow(heat_mat), " x ", ncol(heat_mat), ")")
    return(invisible(NULL))
  }
  with_pdf(file, width = max(8, 0.12 * ncol(heat_mat) + 4), height = 0.18 * nrow(heat_mat) + 4,
    pheatmap::pheatmap(heat_mat, annotation_col = anno, scale = "row", clustering_method = "ward.D2", fontsize_row = 8, fontsize_col = 6,
                       border_color = NA, angle_col = 90, color = pal, main = "Median log2(tpm + 1) per plot_group (row z-score)"))
}
 
## ---------------------------------------------------------------- dot plot
dot_plot <- function(df, group_col, title, feature_labels, min_n = 3, tpm_cut = 0.5, facet_row = NULL, facet_col = "category") {
  dot_df <- df %>%
    pivot_longer(all_of(feature_labels), names_to = "feature", values_to = "expr") %>%
    group_by(across(all_of(c(group_col, facet_row, facet_col))), feature) %>%
    summarise(frac = mean(expr >= log2(tpm_cut + 1)), med = median(expr), n = n(), .groups = "drop") %>%
    filter(n >= min_n)
  if (nrow(dot_df) == 0) return(NULL)
  dot_df$feature <- factor(dot_df$feature, levels = rev(feature_labels))   # first feature at the top
 
  p <- ggplot(dot_df, aes(.data[[group_col]], feature, size = frac, colour = med)) +
    geom_point() +
    scale_colour_gradientn(colours = c("white", "#4575B4", "#D73027"), name = "median of log2(tpm + 1)\nall samples") +
    scale_size(range = c(0.3, 5), limits = c(0, 1), name = paste0("fraction\nTPM >= ", tpm_cut)) +
    labs(x = group_col, y = NULL) +
    my_theme() +
    theme(panel.background = element_rect(fill = "grey92"), panel.grid.major = element_line(colour = "white")) +
    ggtitle(title)
  if (!is.null(facet_row) || !is.null(facet_col)) {
    p <- p + facet_grid(rows = if (!is.null(facet_row)) vars(.data[[facet_row]]),
                        cols = if (!is.null(facet_col)) vars(.data[[facet_col]]),
                        scales = "free_x", space = "free_x")
  }
  p
}
 
plot_dot_pdf <- function(plot_df, feature_labels, file, types, highlight_terms) {
  show_plot <- function(p) if (!is.null(p)) print(p)
  with_pdf(file, width = 25, height = max(5, 0.2 * length(feature_labels) + 3), {
    show_plot(dot_plot(plot_df, "plot_group", "All types combined", feature_labels))
    id_df <- prep_identity(plot_df, highlight_terms)
    if (nrow(id_df) > 0) show_plot(dot_plot(id_df, "identity", "cell-line names", feature_labels, facet_row = "treatment_status", facet_col = NULL))
    for (t in types) {
      sub_df <- plot_df[plot_df$type == t, ]
      if (nrow(sub_df) == 0) next
      show_plot(dot_plot(sub_df, "plot_group", as.character(t), feature_labels))
    }
  })
}
 
## ---------------------------------------------------------------- tau specificity strips
pal_strip <- colorRampPalette(c("#4575B4", "#91BFDB", "#E0F3F8", "#FEE090", "#FC8D59", "#D73027"))(100)  # blue = low, red = high
 
specificity_strip <- function(df, group_col, feature_labels, highlight, min_n = 3, per_feature = FALSE, tile_ratio = 1,
                              title = NULL, zero_cut = 1e-12) {
  grp <- df %>% group_by(.data[[group_col]]) %>%
    summarise(n_samples = n(), across(all_of(feature_labels), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    filter(n_samples >= min_n)
  if (nrow(grp) < 3) return(NULL)
 
  M <- 2^t(as.matrix(grp[, feature_labels])) - 1                  # rows = features, cols = groups, linear TPM
  colnames(M) <- as.character(grp[[group_col]])
  M <- M[rowSums(M) > 0, , drop = FALSE]                          # features never expressed are dropped
  if (nrow(M) == 0) return(NULL)
  tau <- apply(M, 1, function(x) sum(1 - x / max(x)) / (length(x) - 1))
  peak <- apply(M, 1, which.max)
 
  cell <- (M / rowSums(M)) * tau                                  # each feature's tau split across groups by share
  col_order <- names(sort(colSums(cell), decreasing = TRUE))
 
  if (per_feature) {
    rownames(cell) <- sprintf("%s (tau = %.2f)", rownames(cell), tau)
    rord <- order(match(colnames(M)[peak], col_order), -tau)      # staircase
    long <- as.data.frame(as.table(cell[rord, col_order, drop = FALSE]))
    long$Var1 <- factor(long$Var1, levels = rev(rownames(cell)[rord]))
    fill_lab <- "tau x share of\nexpression"
  } else {
    long <- data.frame(Var1 = "all features", Var2 = factor(col_order, levels = col_order), Freq = colSums(cell)[col_order])
    fill_lab <- "sum of\ntau x share"
  }
  long$Freq[abs(long$Freq) < zero_cut] <- NA                      # values below zero_cut become NA
  if (all(is.na(long$Freq))) return(NULL)
  zero_df <- long[is.na(long$Freq), ]                             # cells to draw in grey
  key_lab <- if (zero_cut > 1e-6) paste0("< ", zero_cut) else "0" # legend label for the grey key
 
  # Colour limits = data range; breaks include both ends of the range
  rng <- range(long$Freq, na.rm = TRUE)
  brks <- pretty(rng, n = 5)
  brks <- brks[brks >= rng[1] & brks <= rng[2]]
  gap <- max(diff(rng), 1e-9) * 0.08                              # don't add an end label that would crowd a neighbour
  if (length(brks) == 0 || min(abs(brks - rng[1])) > gap) brks <- c(rng[1], brks)
  if (min(abs(brks - rng[2])) > gap) brks <- c(brks, rng[2])
  brks <- sort(unique(brks))
  brk_labels <- function(x) vapply(x, function(v) format(round(v, 2)), character(1))
 
  hl_col <- if (length(highlight) > 0) ifelse(grepl(paste(highlight, collapse = "|"), col_order, ignore.case = TRUE), "red", "black") else "black"
 
  ggplot(long, aes(Var2, Var1)) +
    geom_tile(aes(fill = Freq)) +
    geom_tile(data = zero_df, fill = "grey90", linewidth = 0) +                       # the grey cells themselves
    geom_point(data = long[1, ], aes(colour = key_lab), alpha = 0, shape = 22) +      # invisible layer: always creates the legend key
    scale_fill_gradientn(colours = pal_strip, name = fill_lab, na.value = "transparent",
                         limits = rng, breaks = brks, labels = brk_labels) +
    scale_colour_manual(values = setNames("grey60", key_lab), name = NULL) +
    guides(fill = guide_colourbar(order = 1),
           colour = guide_legend(order = 2, override.aes = list(shape = 22, fill = "grey90", colour = "grey60", alpha = 1, size = 4))) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = tile_ratio) +
    labs(x = NULL, y = NULL) + ggtitle(title) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, colour = hl_col))
}
 
make_strips <- function(plot_df, feature_labels, types, highlight_terms, per_feature, file, height) {
  with_pdf(file, width = 16, height = height, {
    sp <- function(d, g, ttl) {
      p <- specificity_strip(d, g, feature_labels, highlight = highlight_terms, per_feature = per_feature, title = ttl)
      if (!is.null(p)) print(p)
    }
    sp(plot_df, "plot_group", "All types combined")
    id_df <- prep_identity(plot_df, highlight_terms)
    if (nrow(id_df) > 0) sp(id_df, "identity", "cell-line names")
    for (t in types) sp(plot_df[plot_df$type == t, ], "plot_group", as.character(t))
  })
}

## ---------------------------------------------------------------- expression correlation
# Add grouping factor .grp to df. group_col = "all" pools every sample.
# Existing factor level order is kept (e.g. plot_group follows the highlight order).
cor_group_df <- function(df, group_col = "all") {
  if (identical(group_col, "all")) {
    df$.grp <- factor("all samples")
    return(df)
  }
  if (!group_col %in% colnames(df)) stop("Correlation grouping column not found: ", group_col)
  g <- df[[group_col]]
  df$.grp <- if (is.factor(g)) g else factor(as.character(g))
  df <- df[!is.na(df$.grp) & as.character(df$.grp) != "", , drop = FALSE]
  df$.grp <- droplevels(df$.grp)
  df
}

# Feature x feature correlation matrix for one set of samples (NULL if fewer than min_n samples).
# Constant features have no correlation and stay NA.
cor_matrix <- function(d, features, method = "spearman", min_n = 5) {
  m <- as.matrix(d[, features, drop = FALSE])
  if (nrow(m) < min_n) return(NULL)
  keep <- apply(m, 2, function(v) isTRUE(sd(v, na.rm = TRUE) > 0))
  r <- matrix(NA_real_, length(features), length(features), dimnames = list(features, features))
  if (sum(keep) >= 2) {
    r[keep, keep] <- suppressWarnings(cor(m[, keep, drop = FALSE], method = method, use = "pairwise.complete.obs"))
  }
  r
}

# Long table of all feature pairs per group level: method, group_by, group, feature1, feature2, n, rho, p, padj
# (padj = BH within each group level). Groups with fewer than min_n samples are omitted. NULL if nothing left.
cor_pairs_long <- function(df, features, group_col = "all", method = "spearman", min_n = 5) {
  df <- cor_group_df(df, group_col)
  pairs <- t(combn(features, 2))
  out <- list()
  for (g in levels(df$.grp)) {
    d <- df[df$.grp == g, , drop = FALSE]
    if (nrow(d) < min_n) next
    res <- lapply(seq_len(nrow(pairs)), function(i) {
      x <- d[[pairs[i, 1]]]; y <- d[[pairs[i, 2]]]
      ok <- is.finite(x) & is.finite(y); n <- sum(ok)
      rho <- NA_real_; p <- NA_real_
      if (n >= min_n && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
        rho <- unname(ct$estimate); p <- ct$p.value
      }
      c(n = n, rho = rho, p = p)
    })
    res <- as.data.frame(do.call(rbind, res))
    res$padj <- p.adjust(res$p, method = "BH")
    out[[g]] <- data.frame(method = method, group_by = group_col, group = g,
                           feature1 = pairs[, 1], feature2 = pairs[, 2], res,
                           row.names = NULL, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(NULL)
  out <- do.call(rbind, unname(out))
  rownames(out) <- NULL
  out
}

# Feature x feature heatmap, one page per group level (draws into the open device)
cor_heatmap <- function(df, features, group_col = "all", method = "spearman", min_n = 5,
                        pal = heatmap_palettes$div, max_numbers = 15) {
  df <- cor_group_df(df, group_col)
  for (g in levels(df$.grp)) {
    d <- df[df$.grp == g, , drop = FALSE]
    r <- cor_matrix(d, features, method, min_n)
    if (is.null(r)) next
    usable <- colSums(!is.na(r)) > 0                       # drops constant features
    if (sum(usable) < 2) next
    r <- r[usable, usable, drop = FALSE]
    ord <- hclust(as.dist(1 - r), method = "ward.D2")$order
    r <- r[ord, ord]
    pheatmap::pheatmap(r, color = pal, breaks = seq(-1, 1, length.out = length(pal) + 1),
                       cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey85", border_color = NA,
                       display_numbers = nrow(r) <= max_numbers, number_format = "%.2f",
                       fontsize_row = 8, fontsize_col = 8, angle_col = 90,
                       main = sprintf("%s: %s (n = %d, %s)", group_col, g, nrow(d), method))
  }
}

# Group x group heatmap: how similar are the expression profiles of the groups across the features?
# Median profile per group, centred per feature across groups so baseline expression doesn't dominate.
cor_groups_heatmap <- function(df, features, group_col, method = "spearman", min_n = 5, pal = heatmap_palettes$div) {
  df <- cor_group_df(df, group_col)
  keep <- names(which(table(df$.grp) >= min_n))
  if (length(keep) < 3 || length(features) < 3) return(invisible(NULL))
  prof <- sapply(keep, function(g) apply(df[df$.grp == g, features, drop = FALSE], 2, median, na.rm = TRUE))
  prof <- prof - rowMeans(prof)
  prof <- prof[, apply(prof, 2, sd) > 0, drop = FALSE]    # flat profiles have no correlation
  if (ncol(prof) < 3) return(invisible(NULL))
  r <- suppressWarnings(cor(prof, method = method, use = "pairwise.complete.obs"))
  hc <- hclust(as.dist(1 - r), method = "ward.D2")
  pheatmap::pheatmap(r, color = pal, breaks = seq(-1, 1, length.out = length(pal) + 1),
                     cluster_rows = hc, cluster_cols = hc, border_color = NA, na_col = "grey85",
                     fontsize_row = 6, fontsize_col = 6, angle_col = 90,
                     main = sprintf("%s: expression-profile similarity (%s)", group_col, method))
}

# Dot plot of the top feature pairs (rows) across group levels (columns): fill = rho, size = -log10(padj)
# Input: rows of cor_pairs_long() for one method and one group_by.
cor_dotplot <- function(cor_long, max_pairs = 30) {
  cl <- cor_long[!is.na(cor_long$rho), ]
  if (nrow(cl) == 0) return(NULL)
  cl$pair <- paste(cl$feature1, "~", cl$feature2)
  top <- head(names(sort(tapply(abs(cl$rho), cl$pair, max), decreasing = TRUE)), max_pairs)
  cl <- cl[cl$pair %in% top, ]
  cl$pair <- factor(cl$pair, levels = rev(top))
  cl$group <- factor(cl$group, levels = unique(cor_long$group))
  cl$sig <- pmin(-log10(cl$padj), 10)
  ggplot(cl, aes(group, pair, fill = rho, size = sig)) +
    geom_point(shape = 21, colour = "grey30") +
    scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), limits = c(-1, 1),
                         name = paste0(unique(cl$method), "\nrho")) +
    scale_size_continuous(range = c(0.5, 5), limits = c(0, 10), name = "-log10(padj)") +
    labs(x = unique(cl$group_by), y = NULL) +
    my_theme() +
    theme(panel.background = element_rect(fill = "grey92"), panel.grid.major = element_line(colour = "white")) +
    ggtitle(paste("Top feature-pair correlations by", unique(cl$group_by)))
}

# Scatter plots of the strongest pairs (input: rows of cor_pairs_long() with group_by == "all"); returns a list of ggplots
cor_scatter <- function(df, cor_all, color_col = "type", max_pairs = 12) {
  cl <- cor_all[!is.na(cor_all$rho), ]
  if (nrow(cl) == 0) return(NULL)
  cl <- head(cl[order(-abs(cl$rho)), ], max_pairs)
  lapply(seq_len(nrow(cl)), function(i) {
    d <- df
    d$x <- d[[cl$feature1[i]]]; d$y <- d[[cl$feature2[i]]]
    ggplot(d, aes(x, y, colour = .data[[color_col]])) +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
      scale_colour_viridis_d(name = color_col, drop = TRUE) +
      labs(x = paste(cl$feature1[i], "log2(tpm + 1)"), y = paste(cl$feature2[i], "log2(tpm + 1)")) +
      ggtitle(sprintf("%s rho = %.2f (n = %d)", cl$method[i], cl$rho[i], cl$n[i])) +
      my_theme(axis_title_size = 7, legend_size = 8) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))
  })
}

# Driver: every method x every grouping column.
# Writes, per method m:  <outstem>.cor.<m>.heatmap.pdf   (feature x feature, one page per group level, all group_cols)
#                        <outstem>.cor.<m>.groups.pdf    (group x group profile similarity, one page per group_col)
#                        <outstem>.cor.<m>.dotplot.pdf   (top pairs x groups, one page per group_col)   [pair outputs only if
#                        <outstem>.cor.<m>.scatter.pdf   (strongest pairs over all samples)               n features <= max_pairs_features]
# and one tidy table <outstem>.cor.pairs.tsv with all methods and groupings (ready for the viewer).
plot_correlation <- function(plot_df, feature_labels, outdir, outstem, group_cols = "all",
                             methods = c("spearman", "pearson", "kendall"), min_n = 5,
                             max_pairs_features = 50, max_dot_pairs = 30, max_scatter_pairs = 12,
                             scatter_color = "type") {
  if (length(feature_labels) < 2) { message("Correlation skipped: need >= 2 features"); return(invisible(NULL)) }
  group_cols <- unique(c("all", group_cols))                       # "all" is always computed
  absent <- setdiff(group_cols, c("all", colnames(plot_df)))
  if (length(absent) > 0) message("Correlation: column(s) not found, skipping: ", paste(absent, collapse = ", "))
  group_cols <- setdiff(group_cols, absent)
  n_ok <- vapply(group_cols, function(g) sum(table(cor_group_df(plot_df, g)$.grp) >= min_n), numeric(1))
  if (any(n_ok == 0)) message("Correlation: no group with >= ", min_n, " samples, skipping: ", paste(group_cols[n_ok == 0], collapse = ", "))
  group_cols <- group_cols[n_ok > 0]; n_ok <- n_ok[n_ok > 0]
  if (length(group_cols) == 0) return(invisible(NULL))
  message("Correlation (min ", min_n, " samples/group) by: ", paste0(group_cols, " [", n_ok, " groups]", collapse = ", "))

  prof_cols <- setdiff(group_cols, "all")
  do_pairs  <- length(feature_labels) <= max_pairs_features
  if (!do_pairs) message("More than ", max_pairs_features, " features: skipping pair tables, dot plots and scatter plots")
  hm_size   <- max(8, 0.14 * length(feature_labels) + 4)
  all_pairs <- list()

  for (m in methods) {
    message("  correlation method: ", m)
    stem <- paste0(outdir, "/", outstem, ".cor.", m)
    with_pdf(paste0(stem, ".heatmap.pdf"), hm_size, hm_size,
             for (g in group_cols) cor_heatmap(plot_df, feature_labels, g, m, min_n))
    if (length(prof_cols) > 0) {
      with_pdf(paste0(stem, ".groups.pdf"), 9, 9,
               for (g in prof_cols) cor_groups_heatmap(plot_df, feature_labels, g, m, min_n))
    }
    if (do_pairs) {
      pm <- do.call(rbind, lapply(group_cols, function(g) cor_pairs_long(plot_df, feature_labels, g, m, min_n)))
      if (is.null(pm)) next
      all_pairs[[m]] <- pm
      if (length(prof_cols) > 0) {
        with_pdf(paste0(stem, ".dotplot.pdf"), 14, 0.22 * max_dot_pairs + 3, {
          for (g in prof_cols) {
            p <- cor_dotplot(pm[pm$group_by == g, ], max_dot_pairs)
            if (!is.null(p)) print(p)
          }
        })
      }
      sc <- cor_scatter(plot_df, pm[pm$group_by == "all", ], scatter_color, max_scatter_pairs)
      if (!is.null(sc)) {
        with_pdf(paste0(stem, ".scatter.pdf"), 14, 9, {
          for (s in seq(1, length(sc), by = 6)) print(patchwork::wrap_plots(sc[s:min(s + 5, length(sc))], ncol = 3, guides = "collect"))
        })
      }
    }
  }
  if (length(all_pairs) > 0) {
    write.table(do.call(rbind, unname(all_pairs)), paste0(outdir, "/", outstem, ".cor.pairs.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(all_pairs)
}