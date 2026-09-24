# ============================================================
# FUNCTIONS
# ============================================================
# Written by Isabela Almeida
# Created on May 13, 2026
# Last modified on Sep 24, 2026
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
    write.table(as.data.frame(res) %>% dplyr::select(-leadingEdge), file.path(outdir, paste0(outstem, ".", ctr, ".", db, ".fgsea.tsv")),
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

plot_heatmap <- function(tpmcnt, groups_vec, genesets, qlfo_list,
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
