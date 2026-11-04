args <- commandArgs(trailingOnly = TRUE)
input_counts <- args[which(args == "--inputc") + 1]
input_metadata <- args[which(args == "--inputm") + 1]
input_features <- args[which(args == "--inputf") + 1]
star_dir <- args[which(args == "--stardir") + 1]
highlight <- args[which(args == "--highlight") + 1]
exclude <- args[which(args == "--exclude") + 1]
max_features <- args[which(args == "--maxfeatures") + 1]
outdir <- args[which(args == "--outdir") + 1]
outstem <- args[which(args == "--outstem") + 1]
functions <- args[which(args == "--function") + 1]

## TESTING ONLY: set to TRUE to use the hardcoded values below instead of command-line arguments
testing <- TRUE
if (testing) {
  input_metadata <- "/working/lab_julietF/isabelaA/Project_RNA-Hub_IA/metadata/metadata.tsv"
  #input_counts <- "/working/lab_julietF/isabelaA/Project_HyDRA_IA/rnaseq/grade055_rquant_Bash_04022026141102AEST/RSEM-hydra-rnaseq_tpm-headerfixed.tsv"
  input_counts <- "/working/lab_julietF/isabelaA/Project_BC-OC-lncRNAs_IA/lab_lncRNAs4LW-SC-SE-FC/rnaseq/grade055_rquant_Bash_29052026092303AEST/RSEM-bcoc_tpm.header-fixed.tsv"
  star_dir <- "/working/lab_julietF/isabelaA/Project_RNA-Hub_IA/rnaseq/7kdatarnaseq/grade051_alignment_STAR_09102026200034AEST/"
  #input_features <- "/working/lab_julietF/isabelaA/Project_RNA-Hub_IA/rnaseq/7kdatarnaseq/features.gencode.lncRNAs.transcript.tsv"
  #grep -P "lincRNA|lncRNA" /working/lab_julietF/isabelaA/data/references/TranscriptomeGRCh38rel79.gtf | grep -P "\ttranscript\t" | cut -d'"' -f6,16 | tr '"' '\t' | sort | uniq | while read l; do col1=`echo "${l}" | cut -f1`; col2=`echo "${l}" | cut -f2`; echo -e "${col1}\t${col2}\n${col1}_PAR_Y\t${col2}" ; done > features.gencode.lncRNAs.transcript.tsv
  input_features <- "/working/lab_julietF/isabelaA/Project_RNA-Hub_IA/rnaseq/7kdatarnaseq/features.bcoclab.lncRNAs.transcript.tsv"
  #cut -f1 Project_BC-OC-lncRNAs_IA/lab_lncRNAs4LW-SC-SE-FC/rnaseq/grade055_rquant_Bash_29052026092303AEST/RSEM-bcoc_tpm.header-fixed.tsv | grep -v "ENST" | grep -v "NODE" | grep -v "TRINITY" | while read line ; do echo -e "${line}\t${line}" ; done > /working/lab_julietF/isabelaA/Project_RNA-Hub_IA/rnaseq/7kdatarnaseq/features.bcoclab.lncRNAs.transcript.tsv
  highlight <- "breast, ovar"
  exclude <- "MRC5"    # comma-separated, like highlight; "" for none
  max_features <- 50   # character(0) = default (50)
  #outstem <- "gencode.lncRNAs.transcript"
  outstem <- "bcoclab.lncRNAs.transcript"
  outdir <- "/working/lab_julietF/isabelaA/Project_RNA-Hub_IA/metadata/plots/test/"
  functions <- "/working/lab_julietF/isabelaA/scripts/GRADE2/scripts/ipda_GRADE2_rfunctions.r"
}

print_help <- function() {
  cat("
Written by Isabela Almeida
Created on May 13, 2026
Last modified on Nov 04, 2026
Version: 1.0.0

Description: Plot sample metadata, STAR alignment metrics and TPM expression of a feature list
(violin, heatmap, dot plot and TAU specificity strips) for step 083 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: Rscript ipda_grade2_step083.r [options]

Options:
  --inputc FILE        Input TPM counts TSV file from grade055_rquant_Bash_DATE (needs columns: transcript, gene)
  --inputm FILE        Input metadata TSV file (needs columns: sequencing.id, type, category, biological.system,
                       source, phenotype, phenotype.subtype, phenotype.identity, treatment)
  --inputf FILE        Input features TSV file, no header, 2 columns: feature_id (transcript or gene id), feature_name
  --stardir DIR        [optional] Directory with grade051_alignment_STAR_DATE output; skips STAR metrics if absent
  --highlight STRING   [optional] Comma-separated terms to highlight/order first in phenotype groups (e.g. 'breast,ovar')
  --exclude STRING     [optional] Comma-separated terms; samples with any metadata column matching are dropped (e.g. 'MRC5')
  --maxfeatures INT    [optional] Max number of features for per-feature PDFs and per-feature tau strips (default: 50)
  --outdir DIR         Output directory
  --outstem STEM       Output file stem
  --function FILE      Path to R functions file ipda_grade2_rfunctions.r
  --help               Show this help message

Example:
  Rscript ipda_grade2_step083.r --inputc /path/from/working/dir/to/grade055_rquant_Bash_DATE/RSEM-stem_tpm.tsv --inputm /path/from/working/dir/to/metadata.tsv --inputf /path/from/working/dir/to/features.tsv --stardir /path/from/working/dir/to/grade051_alignment_STAR_DATE --highlight 'breast,ovar' --exclude MRC5 --outdir /path/from/working/dir/to/grade083_plots_R_DATE --outstem stem --function /path/from/working/dir/to/ipda_grade2_rfunctions.r

Pipeline description:

#   000 Index building (0gffcompare, 1Kallisto, 2RSEM, 3STAR, 4Salmon)
#   010 Quality check raw files (0Bedtools, 1FastQC, 2MultiQC)
#   020 Trim reads of adapters (1Trimmomatic)
#   030 Quality check trimmed files (1FastQC, 2MultiQC)
#   040 Pseudo align and quantify reads (1Kallisto, 2BASH count tables)
#   050 Align (1STAR, 2SAMtools, 3NovoSort) and quantify reads (4RSEM, 5BASH count tables)
#   060 PSeudo align and quantify reads at isoform level (1Salmon, 2BASH count tables)
#   070 Differential Expression Analysis (1EdgeR)
#-->083 Plot counts and draw metrics (1pca+heatmap+boxplot, 2metrics, 3id plots)

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
\n")
}

# Show help if requested or no args
if (length(args) == 0 || "--help" %in% args) {
  print_help()
  quit(save = "no")
}

## Check arguments (flag absent or last on the command line gives character(0))
required <- list("--inputc" = input_counts, "--inputm" = input_metadata, "--inputf" = input_features,
                 "--outdir" = outdir, "--outstem" = outstem, "--function" = functions)
missing_args <- names(required)[lengths(required) == 0]
if (length(missing_args) > 0) {
  stop("Missing required argument(s): ", paste(missing_args, collapse = ", "), ". Run with --help for usage.")
}
if (length(star_dir) == 0) star_dir <- NULL
if (length(highlight) == 0) highlight <- ""
if (length(exclude) == 0) exclude <- ""
max_features <- if (length(max_features) == 0) 50 else as.integer(max_features)
for (f in c(input_counts, input_metadata, input_features, functions)) {
  if (!file.exists(f)) stop("File not found: ", f)
}
if (!is.null(star_dir) && !dir.exists(star_dir)) stop("STAR directory not found: ", star_dir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidytext)
  library(patchwork)
  library(data.table)
  library(ggbeeswarm)
  library(pheatmap)
})
needed <- c("paletteer", "scales", "viridis", "RColorBrewer")
not_installed <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed) > 0) stop("Missing R packages: ", paste(not_installed, collapse = ", "))

## Functions
source(functions)

## Settings
highlight_terms <- split_terms(highlight)
exclude_terms <- split_terms(exclude)
type_levels <- c("patient-tissue", "cell-line", "cell-type", "PDX", "CDX", "biofluid")  # plotting order; unknown types are appended
required_meta_cols <- c("sequencing.id", "type", "category", "biological.system", "source", "phenotype",
                        "phenotype.subtype", "phenotype.identity", "treatment")

## Load metadata, features and expression matrix
# Metadata
metadata <- read.delim(input_metadata)
missing_cols <- setdiff(required_meta_cols, colnames(metadata))
if (length(missing_cols) > 0) stop("Metadata is missing column(s): ", paste(missing_cols, collapse = ", "))
metadata <- exclude_samples(metadata, exclude_terms)

# Counts matrix
counts <- fread(input_counts, data.table = FALSE)
if (!all(c("transcript", "gene") %in% colnames(counts))) stop("Counts file needs 'transcript' and 'gene' columns")
counts$TranscriptID <- sub("\\.[0-9]+$", "", counts$transcript)
counts$GeneID <- sub("\\.[0-9]+$", "", counts$gene)

# Features (i.e. transcript/gene ids for plotting)
features <- read.delim(input_features, header = FALSE)
colnames(features)[1:2] <- c("feature_id", "feature_name")
features <- distinct(features, feature_id, .keep_all = TRUE)  # avoid duplicated rows after the join

# Subset by feature list and match samples in counts and metadata
common_samples <- intersect(colnames(counts), metadata$sequencing.id)
if (length(common_samples) == 0) stop("No sample IDs shared between counts columns and metadata$sequencing.id")
message(length(common_samples), " samples shared between counts and metadata")

counts_subset <- counts %>% filter(TranscriptID %in% features$feature_id | GeneID %in% features$feature_id) %>%
  select(transcript, gene, TranscriptID, GeneID, all_of(common_samples))
if (nrow(counts_subset) == 0) stop("None of the features in --inputf were found in the counts file")
metadata <- metadata %>% filter(sequencing.id %in% common_samples)
counts_subset$feature_id <- ifelse(counts_subset$TranscriptID %in% features$feature_id, counts_subset$TranscriptID, counts_subset$GeneID)
counts_subset <- counts_subset %>% left_join(features, by = "feature_id")
counts_subset$feature_label <- paste0(counts_subset$feature_name, " (", counts_subset$feature_id, ")")
feature_labels <- counts_subset$feature_label
message(length(feature_labels), " features matched")

## Process counts data for plotting
# log2(tpm + 1); drop = FALSE keeps this working with a single sample
norm_counts <- log2(counts_subset[, common_samples, drop = FALSE] + 1)
rownames(norm_counts) <- counts_subset$feature_label
tnorm_counts <- as.data.frame(t(norm_counts))
tnorm_counts$sample <- rownames(tnorm_counts)
tnorm_counts <- left_join(tnorm_counts, metadata, by = c("sample" = "sequencing.id"))

plot_df <- add_plot_group(tnorm_counts, highlight_terms, type_levels)
types <- levels(plot_df$type)

## Plot metadata
p_type  <- metadata_barplot_facet(metadata, "type", "category", "biological.system") + ggtitle("Type")
p_pheno <- metadata_barplot_facet(metadata, "phenotype", "category", "biological.system") + ggtitle("Phenotype")
p_subty <- metadata_barplot_facet(metadata, "phenotype.subtype", "category", "biological.system") + ggtitle("Phenotype Subtype")
p_ident <- metadata_barplot_facet(metadata, "phenotype.identity", "category", "biological.system") + ggtitle("Phenotype Identity")
p_source <- metadata_barplot_facet(metadata, "source", "category", "biological.system") + ggtitle("Source")
with_pdf(paste0(outdir, "/", outstem, ".metadata.pdf"), width = 20, height = 5, {
  for (p in list(p_type, p_pheno, p_subty, p_ident, p_source)) print(p)
})

## Plot raw/alignment stats (only when --stardir is given)
if (!is.null(star_dir)) {
  star_qc <- build_star_qc(metadata, star_dir)
  unmapped_df <- star_qc %>%
    select(source, category, unmapped_mismatch_pct, unmapped_short_pct, unmapped_other_pct) %>%
    pivot_longer(starts_with("unmapped"), names_to = "reason", values_to = "pct")
  plots <- list(
    plotscaledin      = histogram_plot(star_qc, "input_reads", wrap = "source~category", title = "Input reads distribution/depth - scaled", scales = "fixed", log10 = TRUE),
    plotfreein        = histogram_plot(star_qc, "input_reads", wrap = "source~category", title = "Input reads distribution/depth - free scales", scales = "free", log10 = TRUE),
    plotscaleduniq    = histogram_plot(star_qc, "unique_reads", wrap = "source~category", title = "Unique reads distribution/depth - scaled", scales = "fixed", log10 = TRUE),
    plotfreeuniq      = histogram_plot(star_qc, "unique_reads", wrap = "source~category", title = "Unique reads distribution/depth - free scales", scales = "free", log10 = TRUE),
    plotscaledmapdist = violin_distribution(star_qc, "source", "unique_map_pct", gridcol = "category", title = "Mapping rate distribution - scaled", scales = "fixed", show_n = FALSE),
    plotfreemapdist   = violin_distribution(star_qc, "source", "unique_map_pct", gridcol = "category", title = "Mapping rate distribution - free scales", scales = "free", show_n = FALSE),
    plotscaledlen     = barplot_distribution(star_qc, "read_length", wrap = "source~category", title = "Read length distribution - scaled", scales = "fixed"),
    plotfreelen       = barplot_distribution(star_qc, "read_length", wrap = "source~category", title = "Read length distribution - free scales", scales = "free"),
    plotscaledunmap   = boxplot_distribution(unmapped_df, "reason", "pct", wrap = "source~category", title = "Unmapped percentages - scaled", scales = "fixed"),
    plotfreeunmap     = boxplot_distribution(unmapped_df, "reason", "pct", wrap = "source~category", title = "Unmapped percentages - free scales", scales = "free"),
    plotscaledchim    = histogram_plot(star_qc, "chimeric_pct", wrap = "source~category", title = "Chimeric reads - scaled", scales = "fixed", log10 = FALSE),
    plotfreechim      = histogram_plot(star_qc, "chimeric_pct", wrap = "source~category", title = "Chimeric reads - free scales", scales = "free", log10 = FALSE),
    plotmapdepth      = scatter_plot(star_qc, "input_reads", "unique_map_pct", gridcol = "category", title = "Mapping rate vs read depth", log10_x = TRUE, shape = "source")
  )
  plot_heights <- c(plotscaledmapdist = 5, plotfreemapdist = 5, plotmapdepth = 5)  # every other plot: 10
  heights <- setNames(as.list(ifelse(names(plots) %in% names(plot_heights), plot_heights[names(plots)], 10)), names(plots))
  save_plot_list_pdf(plots, heights, paste0(outdir, "/", outstem, ".star.alignment.pdf"))
} else {
  message("--stardir not provided: skipping STAR alignment plots")
}

## Plot expression matrix - features
# Violins of mean expression across all features
plot_mean_pdf(plot_df, feature_labels, paste0(outdir, "/", outstem, ".violin.mean.pdf"), types, highlight_terms)

# One PDF per feature, only for small feature sets to avoid too many files
if (length(feature_labels) <= max_features) {
  file_stems <- sub(" \\(.*$", "", feature_labels)          # "MALAT1 (ENSG...)" -> "MALAT1"
  file_stems <- gsub("[^A-Za-z0-9._-]+", "_", file_stems)   # make filesystem-safe
  file_stems <- make.unique(file_stems, sep = "_")          # avoid overwriting if names repeat
  for (i in seq_along(feature_labels)) {
    plot_feature_pdf(feature_labels[i], file_stems[i], plot_df, outdir, outstem, types, highlight_terms)
  }
} else {
  message("More than ", max_features, " features (", length(feature_labels), "): skipping per-feature PDFs")
}

# Heatmap (median per plot_group)
plot_heatmap(plot_df, feature_labels, paste0(outdir, "/", outstem, ".heatmap.pdf"), pal = heatmap_palettes$div_blue)

# Dot plot of fraction of samples expressing vs median expression
plot_dot_pdf(plot_df, feature_labels, paste0(outdir, "/", outstem, ".dotplot.pdf"), types, highlight_terms)

# Specificity score - TAU
make_strips(plot_df, feature_labels, types, highlight_terms, per_feature = FALSE,
            file = paste0(outdir, "/", outstem, ".tau.combined.pdf"), height = 4)  # Aggregate strips (all features pooled)
if (length(feature_labels) <= max_features) {                                      # One row per feature, only for small feature sets
  make_strips(plot_df, feature_labels, types, highlight_terms, per_feature = TRUE,
              file = paste0(outdir, "/", outstem, ".tau.features.pdf"), height = 0.2 * length(feature_labels) + 4)
}

message("Step 083 finished. Outputs written to: ", outdir)