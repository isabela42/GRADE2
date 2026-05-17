args <- commandArgs(trailingOnly = TRUE)
input_counts <- args[which(args == "--inputc") + 1]
input_metadata <- args[which(args == "--inputm") + 1]
plot_condition <- args[which(args == "--condition") + 1]
plot_sec_condition <- args[which(args == "--seccondition") + 1]
outdir <- args[which(args == "--outdir") + 1]
outstem <- args[which(args == "--outstem") + 1]
functions <- args[which(args == "--function") + 1]

print_help <- function() {
  cat("
Written by Isabela Almeida with input from Larissa Cassiano
Created on May 13, 2026
Last modified on May 17, 2026
Version: 1.0.0

Description: Write and submit PBS jobs for step 055 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: Rscript ipda_camels_step083.r [options]

Options:
  --inputc FILE          Input TPM counts TSV file from grade055_rquant_Bash_DATE
  --inputm FILE          Input metadata TSV file with condition and seccondition columns per sample
  --condition STRING     Input per_target summary TSV file from camels034_combined-summary_BASH-R_DATE
  --seccondition STRING  Input per_target summary TSV file from camels034_combined-summary_BASH-R_DATE
  --outdir DIR           Output directory
  --outstem STEM         Output file stem
  --function FILE        Path to R functions file ipda_grade2_rfunctions.r
  --help                 Show this help message

Example:
  Rscript ipda_grade2_step083.r --inputc /path/from/working/dir/to/grade055_rquant_Bash_DATE/RSEM-stem_tpm.tsv --inputm /path/from/working/dir/to/metadata.tsv --condition string1-from-metadata --seccondition string2-from-metadata --outdir /path/from/working/dir/to/grade083_gene-plots_R_DATE --outstem stem --function /path/from/working/dir/to/ipda_grade2_rfunctions.r

Pipeline description:

#   000 Index building (0gffcompare, 1Kallisto, 2RSEM, 3STAR, 4Salmon)
#   010 Quality check raw files (0Bedtools, 1FastQC, 2MultiQC)
#   020 Trim reads of adapters (1Trimmomatic)
#   030 Quality check trimmed files (1FastQC, 2MultiQC)
#   040 Pseudo align and quantify reads (1Kallisto, 2BASH count tables)
#   050 Align (1STAR, 2SAMtools, 3NovoSort) and quantify reads (4RSEM, 5BASH count tables)
#   060 PSeudo align and quantify reads at isoform level (1Salmon, 2BASH count tables)
#   070 Differential Expression Analysis (1EdgeR)
#-->080 Plot counts and draw metrics (1pca+heatmap+boxplot, 2metrics, 3id plots)

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
\n")
}

# Show help if requested or no args
if (length(args) == 0 || "--help" %in% args) {
  print_help()
  quit(save = "no")
}

## Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

## Set input/output paths
out_boxplot <- file.path(outdir, paste0(outstem, ".", plot_condition, ".boxplot.pdf"))
out_beeswarm <- file.path(outdir, paste0(outstem, ".", plot_condition, ".beeswarm.pdf"))
out_heatmap <- file.path(outdir, paste0(outstem, ".", plot_condition, ".heatmap.pdf"))

## Import data/metadata
counts <- fread(file.path(input_counts), data.table=FALSE)  # Use fread for speed if data.table available
is_transcript <- outstem %in% counts$transcript
is_gene <- outstem %in% counts$gene

if (!is_transcript && !is_gene) {
  stop(
    paste0(
      "'", outstem,
      "' was not found in either transcript or gene columns."
    )
  )
}

if (is_transcript) {
  counts <- counts %>%
    filter(transcript == outstem)
} else if (is_gene) {
  counts <- counts %>%
    filter(gene == outstem)
  sample_cols <- setdiff(
    colnames(counts),
    c("transcript", "gene")
  )
  summed_counts <- colSums(
    counts[, sample_cols, drop = FALSE],
    na.rm = TRUE
  )
  counts <- data.frame(
    transcript = outstem,
    gene = outstem,
    t(summed_counts),
    check.names = FALSE
  )
}
counts <- counts %>% select(-2) # gene column
#counts2 <- counts %>% filter(transcript == outstem) # with the conditioning statement above user can provide either transcript or gene ID
metadata <- read.delim(file.path(input_metadata), sep="\t", header=TRUE, row.names=1)
metadata[[plot_condition]] <- factor(metadata[[plot_condition]]) # Ensure condition is factor

## Process data for plotting
# Subset - match samples in counts and metadata
common_samples <- intersect(colnames(counts)[-1], rownames(metadata)) # Skip gene column - usefull for subsetting data
counts_subset <- counts[, c("transcript", common_samples), drop = FALSE] # Only shared samples
metadata <- metadata[common_samples, , drop = FALSE]

# Calculate log2(tpm + 1) for plotting
norm_counts <- log2(counts_subset[, -1] + 1)
#norm_counts <- counts_subset[, -1]
rownames(norm_counts) <- counts_subset$transcript
tnorm_counts <- as.data.frame(t(norm_counts))
tnorm_counts$samples <- rownames(tnorm_counts)
metadata$samples <- rownames(metadata)
tnorm_counts <- left_join(tnorm_counts, metadata, by = c("samples" = "samples"))
plot_heat <- tnorm_counts %>%
    group_by(.data[[plot_condition]], .data[[plot_sec_condition]]) %>%
    summarise(expr = mean(.data[[outstem]], na.rm = TRUE), .groups = "drop")
plot_heat <- plot_heat %>%
  complete(
    !!sym(plot_sec_condition),
    !!sym(plot_condition)
  )
## Source functions
source(functions)

## Boxplot
boxplot <- box_plot(tnorm_counts, plot_condition, outstem)
ggsave(file.path(out_boxplot), plot = boxplot, width = 20, height = 5, dpi = 100)

## Beeswarm
beeswarm <- beeswarm_plot(tnorm_counts, plot_condition, outstem)
ggsave(file.path(out_beeswarm), plot = beeswarm, width = 20, height = 5, dpi = 100)

## Heatmap
heatmap <- heatmap_plot(plot_heat)
ggsave(file.path(out_heatmap), plot = heatmap, width = 20, height = 20, dpi = 100)
