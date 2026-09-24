args <- commandArgs(trailingOnly = TRUE)
input_anno <- args[which(args == "--inputanno") + 1]
input_expcnt <- args[which(args == "--inputexp") + 1]
input_tpmcnt <- args[which(args == "--inputtpm") + 1]
level <- args[which(args == "--level") + 1]
groups <- args[which(args == "--groups") + 1]
colnums <- args[which(args == "--colnums") + 1]
contrasts <- args[which(args == "--contrasts") + 1]
merge <- args[which(args == "--merge") + 1]
mergecontrasts <- args[which(args == "--mergecontrasts") + 1]
genesets <- args[which(args == "--genesets") + 1]
plotlabel <- args[which(args == "--plotlabel") + 1]
input_gmt <- args[which(args == "--inputgmt") + 1]
outdir <- args[which(args == "--outdir") + 1]
outstem <- args[which(args == "--outstem") + 1]
functions <- args[which(args == "--function") + 1]

print_help <- function() {
  cat("
Written by Isabela Almeida
Based on
  - Larissa Cassiano's script
Created on Sep 11, 2026
Last modified on Sep 24, 2026
Version: 1.0.0

Description: Write and submit PBS jobs for step 071 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: Rscript ipda_grade2_step071.r [options]

Options:
  --inputanno FILE         Input annotation file with GENE_ID|TRANSCRIPT_ID|GENE_SYMBOL|TRANSCRIPT_SYMBOL where IDs match those of inputexp and inputtpm (header false)
  --inputexp FILE          Input expected counts TSV file from grade055_rquantBash_DATE
  --inputtpm FILE          Input TPM counts TSV file from grade055_rquant_Bash_DATE
  --level STRING           gene or transcript (which column to use from counts file)
  --groups STRING          Comma separated group labels based on order of EXP/TPM counts file, e.g.: \"CTRL,ASO1,ASO2,CTRL,ASO1,ASO2,CTRL,ASO1,ASO2\"
  --colnums STRING         Column numbers from inputexp and inputtpm to use. Always include 1,2 (transcript,gene) followed by the other column numbers.
  --contrasts STRING       Pairwise comma separated raw contrast groups for differential expression analysis, e.g.: \"ASO1_vs_CTRL,ASO2_vs_CTRL\"
  --merge STRING           Merging rule for group labels, e.g.: \"ASO1,ASO2=ASO\"
  --mergecontrasts STRING  Pairwise comma separated merged contrast groups for differential expression analysis, e.g.: \"ASO_vs_CTRL\"
  --genesets STRING        Highlight inbuilt gene sets list on plots (use \"NA\" for none, options: HR,FA,BER,NHEJ,Cell_Cycle,top10DEG)
  --plotlabel STRING       Label for plots
  --inputgmt FILE          Input GSEA MSigDB collections folder, e.g. /path/from/working/dir/to/msigdb_collections/human - download from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
  --outdir DIR             Output directory
  --outstem STEM           Output file stem
  --function FILE          Path to R functions file ipda_grade2_rfunctions.r
  --help                   Show this help message

Example:
  Rscript ipda_grade2_step071.r --inputanno /path/from/working/dir/to/annotation.tsv --inputexp /path/from/working/dir/to/grade055_rquant_Bash_DATE --inputtpm /path/from/working/dir/to/grade055_rquant_Bash_DATE --level transcript --groups \"CTRL,ASO1,ASO2,CTRL,ASO1,ASO2,CTRL,ASO1,ASO2\" --colnumns \"1,2,12,13,14,15,16,17,18,19,20\" --contrasts \"ASO1_vs_CTRL,ASO2_vs_CTRL\" --merge \"ASO1,ASO2=ASO\" --mergecontrasts \"ASO_vs_CTRL\" --genesets \"HR,FA,BER,NHEJ,Cell_Cycle,top10DEG\" --plotlabel \"top10DEG\" --inputgmt /path/from/working/dir/to/msigdb_collections/human --outdir /path/from/working/dir/to/grade071_DE_R_DATE --outstem stem --function /path/from/working/dir/to/ipda_grade2_rfunctions.r

Pipeline description:

#   000 Index building (0gffcompare, 1Kallisto, 2RSEM, 3STAR, 4Salmon)
#   010 Quality check raw files (0Bedtools, 1FastQC, 2MultiQC)
#   020 Trim reads of adapters (1Trimmomatic)
#   030 Quality check trimmed files (1FastQC, 2MultiQC)
#   040 Pseudo align and quantify reads (1Kallisto, 2BASH count tables)
#   050 Align (1STAR, 2SAMtools, 3NovoSort) and quantify reads (4RSEM, 5BASH count tables)
#   060 PSeudo align and quantify reads at isoform level (1Salmon, 2BASH count tables)
#-->070 Differential Expression Analysis (1EdgeR)
#   080 Plot counts and draw metrics (1pca+heatmap+boxplot, 2metrics, 3id plots)

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
\n")
}

# Show help if requested or no args
if (length(args) == 0 || "--help" %in% args) {
  print_help()
  quit(save = "no")
}

## Load libraries
library(optparse)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(pheatmap)
library(matrixStats)
#library(org.Hs.eg.db)
#library(AnnotationDbi)
library(fgsea)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

## Groups, sets, labels and thresholds
groups_raw      <- strsplit(groups, ",")[[1]]
contrasts_list  <- strsplit(contrasts, ",")[[1]]
merge_groups    <- if (is.na(merge) || merge == "NA") NULL else merge
merge_contrasts <- if (is.na(mergecontrasts) || mergecontrasts == "NA") NULL else strsplit(mergecontrasts, ",")[[1]]
genesets        <- if (is.na(genesets) || genesets == "NA") NULL else strsplit(genesets, ",")[[1]]
plotlabel       <- if (is.na(plotlabel) || plotlabel == "NA") "top5DEG" else plotlabel
min_cpm         <- 3
min_samples     <- 2
fdr_cut         <- 0.5
lfc_cut         <- 1

## Import data
keep_cols <- as.numeric(unlist(strsplit(colnums, split = ",")))
header <- scan(input_expcnt, what = character(), sep = "\t", nlines = 1, quiet = TRUE)
col_types <- rep("NULL", length(header))
col_types[keep_cols] <- NA
expcnt <- read.table(input_expcnt, header=TRUE, sep="\t", colClasses = col_types)
tpmcnt      <- read.table(input_tpmcnt, header=TRUE, sep="\t", colClasses = col_types)
annotation  <- read.table(input_anno, header=FALSE, sep="\t", stringsAsFactors = FALSE)
ngenes      <- length(unique(expcnt[[level]]))

if (level == "transcript") {
    # exp counts
    rownames(expcnt) <- expcnt[[1]]          # transcript_id
    expcnt           <- expcnt[, -c(1,2)]    # remove transcript + gene cols
    
    # tpm counts
    rownames(tpmcnt) <- tpmcnt[[1]]
    tpmcnt           <- tpmcnt[, -c(1,2)]
} else if (level == "gene") {
    # exp counts
    id_col <- colnames(expcnt)[1]
    dupes  <- duplicated(expcnt[[id_col]])
    if (any(dupes))
    expcnt           <- expcnt[!dupes, ]
    rownames(expcnt) <- expcnt[[id_col]]
    expcnt[[id_col]] <- NULL

    # tpm counts
    id_col <- colnames(tpmcnt)[1]
    dupes  <- duplicated(tpmcnt[[id_col]])
    if (any(dupes))
    tpmcnt           <- tpmcnt[!dupes, ]
    rownames(tpmcnt) <- tpmcnt[[id_col]]
    tpmcnt[[id_col]] <- NULL
}
counts_raw <- expcnt

## Merge groups
groups_merged <- groups_raw
if (!is.null(merge_groups)) {
  for (mg in strsplit(merge_groups, ";")[[1]]) {
    parts <- strsplit(mg, "=")[[1]]
    from  <- strsplit(parts[1], ",")[[1]]
    to    <- parts[2]
    groups_merged[groups_merged %in% from] <- to
  }
}

## Annotation
feat_ids <- rownames(counts_raw)
if (level == "transcript") {
    anno <- annotation[match(feat_ids, annotation[, 2]), c(1, 2, 3, 4)]
    colnames(anno) <- c("GENE", "TRANSCRIPT", "SYMBOLGENE", "SYMBOLTRANSCRIPT")
} else {
    anno <- annotation[match(feat_ids, annotation[, 1]), c(1, 3)]
    colnames(anno) <- c("GENE", "SYMBOL")
    anno <- anno[!duplicated(anno$GENE), ]
}

feat_idstpm <- rownames(tpmcnt)
if (level == "transcript") {
    annotpm <- annotation[match(feat_idstpm, annotation[, 2]), c(1, 2, 3, 4)]
    colnames(annotpm) <- c("GENE", "TRANSCRIPT", "SYMBOLGENE", "SYMBOLTRANSCRIPT")
} else {
    annotpm <- annotation[match(feat_idstpm, annotation[, 1]), c(1, 3)]
    colnames(annotpm) <- c("GENE", "SYMBOL")
    annotpm <- annotpm[!duplicated(annotpm$GENE), ]
}

# NOTE: you have the option of using the annotation section below to use AnnotationDbi/org.Hs.eg.db instead. Please note that:
# - this is only possible if you using reference only genes/transcripts
# - if your GTF does not match the same reference version of AnnotationDbi, there may be some IDs which will be assigned to an incorrect gene name
# if (level == "transcript") {
#   feat_ids <- sub("\\.[0-9]+$", "", rownames(counts_raw))
#   anno <- tryCatch(
#     AnnotationDbi::select(org.Hs.eg.db,
#       keys=unique(feat_ids), columns=c("SYMBOL","ENTREZID","ENSEMBL"),
#       keytype="ENSEMBLTRANS"),
#     error=function(e) data.frame(ENSEMBLTRANS=feat_ids, SYMBOL=NA_character_,
#                                   ENTREZID=NA_character_, ENSEMBL=NA_character_)
#   )
#   colnames(anno)[1] <- "ENSEMBL"  # rename for downstream compatibility
# } else {
#   feat_ids <- sub("\\.[0-9]+$", "", rownames(counts_raw))
#   anno <- tryCatch(
#     AnnotationDbi::select(org.Hs.eg.db,
#       keys=unique(feat_ids), columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL"),
#     error=function(e) data.frame(ENSEMBL=feat_ids, SYMBOL=NA_character_,
#                                   ENTREZID=NA_character_)
#   )
# }
# anno <- anno[!duplicated(anno$ENSEMBL), ]


## DNA repair / Cell Cycle gene lists
dna_repair_genes <- list(
    HR         = c("BRCA1","BRCA2","RAD51","RAD51C","RAD51D","RAD54L","PALB2","BARD1",
                    "BRIP1","MRE11A","NBN","RAD50","RBBP8","ATM","ATR","CHEK1",
                    "TOPBP1","CLSPN","RPA1"),
    FA         = c("FANCA","FANCD2","FANCI","FANCM","SLX4","UBE2T"),
    BER        = c("PARP1","XRCC1","PNKP","LIG3","POLB","APEX1"),
    NHEJ       = c("PRKDC","XRCC5","XRCC6","LIG4","NHEJ1","TP53BP1"),
    Cell_Cycle = c("MKI67","PCNA","MCM2","MCM4","MCM7","CCND1","CCNE1","CCNA2",
                    "CCNB1","CDK1","CDK2","CDC25A","WEE1","FOXM1","PLK1","AURKA",
                    "RRM2","E2F1","E2F2","MYBL2"),
    inhouse = c("NODE_518103_length_493_cov_38.995305_g333637_i0", "chr17:32205199-32206771", "BRRIAR", "CUPID1", "CUPID2", "KILLR", "RHOT1", "ATAD5")
)
pathway_colors <- c(HR="#E69F00", FA="#56B4E9", BER="#009E73",
                    NHEJ="#CC79A7", Cell_Cycle="#D55E00", inhouse="blue")

## Source functions
source(functions)
# Note: Use annotate_qlfo_ensemblr(qlft, anno) for AnnotationDbi/org.Hs.eg.db-based annotation

## Run contrasts
all_qlfo <- list()

for (ctr in contrasts_list) {
    res  <- run_edger(counts_raw, groups_raw, ctr)
    qlfo <- annotate_qlfo(res$qlft, anno, level)
    stem <- paste0(outstem, "_", level, "_", ctr)
    write.table(qlfo, file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".full-table.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(qlfo[qlfo$diffexpressed != "NO", ], file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".DE-all.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(qlfo[qlfo$diffexpressed == "UP", ], file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".DE-up.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(qlfo[qlfo$diffexpressed == "DOWN", ], file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".DE-down.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    plot_pca(res$y, groups_raw[groups_raw %in% strsplit(ctr,"_vs_")[[1]]], ctr, outstem, outdir, level)
    plot_volcano(qlfo, ctr, outstem, outdir, volcano_label, level)
    run_fgsea(qlfo, ctr, outstem, input_gmt, outdir, level)
    parts <- strsplit(ctr, "_vs_")[[1]]
    g_subs <- groups_raw[groups_raw %in% parts]
    keep <- groups_raw %in% parts
    tpm_subs <- tpmcnt[, keep, drop = FALSE]
    plot_heatmap(tpm_subs, g_subs, genesets, list(qlfo), outstem, level, outdir, annotpm, ctr)
    all_qlfo[[ctr]] <- qlfo
}

for (ctr in merge_contrasts) {
    res  <- run_edger(counts_raw, groups_merged, ctr)
    qlfo <- annotate_qlfo(res$qlft, anno, level)
    stem <- paste0(outstem, "_", level, "_", ctr, "_merged")
    write.table(qlfo, file.path(outdir, paste0(outstem, ".", level, ".", ctr,".full-table.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(qlfo[qlfo$diffexpressed != "NO", ], file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".DE-all.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(qlfo[qlfo$diffexpressed == "UP", ], file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".DE-up.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(qlfo[qlfo$diffexpressed == "DOWN", ], file.path(outdir, paste0(outstem, ".", level, ".", ctr, ".DE-down.tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    parts  <- strsplit(ctr,"_vs_")[[1]]
    g_subs <- groups_merged[groups_merged %in% parts]
    plot_pca(res$y, g_subs, ctr, outstem, outdir, level)
    plot_volcano(qlfo, ctr, outstem, outdir, volcano_label, level)
    run_fgsea(qlfo, ctr, outstem, input_gmt, outdir, level)
    all_qlfo[[paste0(ctr,"_merged")]] <- qlfo
}

plot_heatmap(tpmcnt, groups_raw, genesets, all_qlfo, outstem, level, outdir, annotpm, ctr)


# ── Save environment ───────────────────────────────────────────────────────────
ts       <- format(Sys.time(), "%d%m%Y%H%M%S")
env_file <- file.path(outdir, paste0("REnvironment_", outstem, ".", level, ".", ts, ".RData"))
save.image(env_file)
tryCatch(savehistory(file.path(outdir,
  paste0("Rhistory_", outstem, ".", level, ".", ts, ".Rhistory"))),
  error=function(e) NULL)

cat("\n## Done!", format(Sys.time()), "\n")