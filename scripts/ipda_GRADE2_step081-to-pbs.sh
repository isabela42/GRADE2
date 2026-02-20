#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Larissa Cassiano's PCA R script
Created on 20 Feb, 2026
Last modified on 20 Feb, 2026
Version: ${version}

Description: Write and submit PBS jobs for step 081 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_081-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources baseline: -m 15 -c 1 -w "00:30:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            plot-name
                            Plots will be saved on grade081_pca_R_DATE/plot-name.pdf

                            Col2:
                            /path/from/working/dir/to/grade*_DATE/*-stem_*.tsv
                            e.g. /path/from/working/dir/to/gradegrade055_rquant_Bash_DATE/RSEM-stem_tpm.tsv

                            Col3:
                            /path/from/working/dir/to/samples-metadata.tsv
                            Note: this file must contain all samples to plot in col1, i.e. 
                            not necessairly all samples in Col2 file followed by a series
                            of columns with conditions to plot - could be just one, could
                            be several. Please note condition to plot must be specified in Col4

                            Col4:
                            condition
                            String of condition to plot, e.g. phenotype

                            Col5:
                            INT
                            Top variable transcripts/genes to select for plot. Recommended a range between 1000-5000

                            Col6:
                            SHAPES
                            For example: 1, 2, 0, 5, 6, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14
                            See detailed explanation below.

-p <PBS stem>               Stem for PBS file names
-e <email>                  Email for PBS job
-m <INT>                    Memory INT required for PBS job in GB
-c <INT>                    Number of CPUS required for PBS job
-w <HH:MM:SS>               Clock walltime required for PBS job

## Output:

logfile.txt                 Logfile with commands executed and date
PBS files                   PBS files created

Pipeline description:

#   000 Index building (0gffcompare, 1Kallisto, 2RSEM, 3STAR, 4Salmon)
#   010 Quality check raw files (0Bedtools, 1FastQC, 2MultiQC)
#   020 Trim reads of adapters (1Trimmomatic)
#   030 Quality check trimmed files (1FastQC, 2MultiQC)
#   040 Pseudo align and quantify reads (1Kallisto, 2BASH count tables)
#-->050 Align (1STAR, 2SAMtools, 3NovoSort) and quantify reads (4RSEM, 5BASH count tables)
#   060 PSeudo align and quantify reads at isoform level (1Salmon, 2BASH count tables)
#   070 Differential Expression Analysis (1EdgeR)
#   081 Plot counts (1pca, 2heatmap)

Example of variations of shapes to use - col6 with one line info, comm+space sepparated:
16, 17, 15, 18, 1, 2, 0, 5, 6, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, \"*\"
1, 2, 0, 5, 6, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14

The is a limitation on the number of shapes to use. These are 80 examples of different shapes:
0, 1, 2, 5, 6,                                       # basic open shapes - 0 square; 1 circle; 2 triangle up; 
                                                     # 5 diamond; 6 triangle down
3, 4, 7, 8, 9, 10, 11, 12, 13, 14,                   # line symbols - 3 +; 4 ×; 7 square×; 8 star;
                                                     # 9 diamond+; 10 circle+; 11 triangle up/down;
                                                     # 12 square+; 13 circle×; 14 filled triangle/square
15, 16, 17, 18, 19, 20,                              # filled solid shapes - 15 square; 16 circle; 17 triangle up;
                                                     # 18 diamond; 19 solid small circle; 20 bullet
21, 22, 23, 24, 25,                                  # filled with outline - 21 circle; 22 square; 23 diamond; 
                                                     # 24 triangle up; 25 triangle down
\"+\", \"-\", \"*\", \".\", \"o\", \"x\", \"#\", \"%\", \"&\", \"=\",   # ASCII symbols - plus, minus, star, dot,
                                                     # circle, cross, hash, percent, ampersand, equals
\"@\", \"[\", \"]\", \"{\", \"}\"                              # extra ASCII symbols

\"A\", \"B\", \"C\", \"D\", \"E\", \"F\", \"G\", \"H\", \"I\", \"J\",    # uppercase letters A-J
\"K\", \"L\", \"M\", \"N\", \"O\", \"P\", \"Q\", \"R\", \"S\", \"T\",    # uppercase letters K-T
\"U\", \"V\", \"W\", \"X\", \"Y\", \"Z\",                        # uppercase letters U-Z
\"0\", \"1\", \"2\", \"3\", \"4\", \"5\", \"6\", \"7\", \"8\", \"9\"     # digits
\"\u25B6\", \"\u25C0\", \"\u27B5\"                         # Unicode - right triangle, left triangle, arrow

Please contact Isabela Almeida at mb.isabela42@gmail.com if you encounter any problems.
"
}

## Display help message
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

#  ____       _   _   _                 
# / ___|  ___| |_| |_(_)_ __   __ _ ___ 
# \___ \ / _ \ __| __| | '_ \ / _` / __|
#  ___) |  __/ |_| |_| | | | | (_| \__ \
# |____/ \___|\__|\__|_|_| |_|\__, |___/
#                             |___/     

## Save execution ID
pid=`echo $$` #$BASHPID

## User name within your cluster environment
user=`whoami`

## Exit when any command fails
set -e

#................................................
#  Input parameters from command line
#................................................

## Get parameters from command line flags
while getopts "i:p:e:m:c:w:h:v" flag
do
    case "${flag}" in
        i) input="${OPTARG}";;       # Input files including path
        p) pbs_stem="${OPTARG}";;    # Stem for PBS file names
        e) email="${OPTARG}";;       # Email for PBS job
        m) mem="${OPTARG}";;         # Memory required for PBS job
        c) ncpus="${OPTARG}";;       # Number of CPUS required for PBS job
        w) walltime="${OPTARG}";;    # Clock walltime required for PBS job
        h) Help ; exit;;             # Print Help and exit
        v) echo "${version}"; exit;; # Print version and exit
        ?) echo script usage: bash ipda_GRADE2_081-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
           exit;;
    esac
done

#................................................
#  Set Logfile stem
#................................................

## Set Logfile stem
# it contains all the executed commands with date/time;
# the output files general metrics (such as size);
# and memory/CPU usage for all executions
thislogdate=$(date +'%d%m%Y%H%M%S%Z')
human_thislogdate=`date`
logfile=logfile_ipda_GRADE2_081-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# RStudio 4.5.0
module_R=R/4.5.0

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
outpath_GRADE2081_R="grade081_pca_R_${thislogdate}"

## Create output directories
mkdir -p ${outpath_GRADE2081_R}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_GRADE2_081-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input folder of files:       ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${input}"
echo "## logfile will be saved as:    ${logfile}"
echo

# ____  _             _   _                   
#/ ___|| |_ __ _ _ __| |_(_)_ __   __ _       
#\___ \| __/ _` | '__| __| | '_ \ / _` |      
# ___) | || (_| | |  | |_| | | | | (_| |_ _ _ 
#|____/ \__\__,_|_|   \__|_|_| |_|\__, (_|_|_)
#                                 |___/  

#................................................
#  Print Execution info to logfile
#................................................

exec &> "${logfile}"

date
echo "## Executing bash ipda_GRADE2_081-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input folder of files:       ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${input}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read plot; do echo "#!/bin/sh" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "##########################################################################" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Script:  ${pbs_stem}_${plot}_${thislogdate}.pbs" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Version: v01" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Email:   ${email}" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "##########################################################################" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read plot; do echo "#PBS -N ${pbs_stem}_${plot}_${thislogdate}" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#PBS -r n" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#PBS -m abe" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#PBS -M ${email}" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Set main working directory" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "## Change to main directory" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo 'echo "## Load tools"' >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "module load ${module_R}" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done

## Write R file
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Write R file" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Writing file to >> ${pbs_stem}_${plot}_${thislogdate}.r" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# PCA and Heatmap/Hclust plots for RNA-seq TPM data" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Load required libraries" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(dplyr)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(data.table)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(ggplot2)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(pheatmap)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(RColorBrewer)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(ggrepel)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(cluster)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(factoextra)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(dendextend)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(DESeq2)  # For rlog/vst (if counts); skip if already TPM log2" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(matrixStats)  # For rowVars" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "library(viridis)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Set input files directory" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "dir <- \"${PBS_O_WORKDIR}\"" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# 1. LOAD DATA" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Read only gene names + top columns first, then subset" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do counts=`grep "${plot}" ${input} | cut -f2 | sort | uniq`; echo "counts <- read.delim(file.path(\"${counts}\"), sep=\"\t\", header=TRUE, nrows=5)  # Check structure" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
#counts <- read.delim(file.path(dir, "grade055_rquant_Bash_04022026141102AEST/RSEM-hydra-rnaseq_tpm-headerfixed.tsv"), sep="\t", header=TRUE, nrows=5)  # Check structure
cut -f1 ${input} | sort | uniq | while read plot; do counts=`grep "${plot}" ${input} | cut -f2 | sort | uniq`; echo "counts <- fread(file.path(\"${counts}\"), data.table=FALSE)  # Use fread for speed if data.table available" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
#counts <- fread(file.path(dir, "grade055_rquant_Bash_04022026141102AEST/RSEM-hydra-rnaseq_tpm-headerfixed.tsv"), data.table=FALSE)  # Use fread for speed if data.table available
cut -f1 ${input} | sort | uniq | while read plot; do echo "counts <- counts %>% select(-2) # gene column" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
#counts_enst <- counts_full %>% filter(transcript %>% grepl("^ENST", .)) #counts_hydra <- counts_full %>% filter(transcript %>% grepl("^NODE", .))
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Load sample metadata" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do metadata=`grep "${plot}" ${input} | cut -f3 | sort | uniq`; echo "coldata <- read.delim(file.path(\"${metadata}\"), sep=\"\t\", header=TRUE, row.names=1)  # rownames = sample names" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do condition=`grep "${plot}" ${input} | cut -f4 | sort | uniq`; echo "coldata\$${condition} <- factor(coldata\$${condition})  # Ensure condition is factor" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Match colnames(counts) with rownames(coldata)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "common_samples <- intersect(colnames(counts)[-1], rownames(coldata))  # Skip gene column - usefull for subsetting data" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "counts_subset <- counts[, c(\"transcript\", common_samples)]  # Only shared samples" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "coldata <- coldata[common_samples, ]" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "cat(\"Samples kept:\", length(common_samples), \"\n\")" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# 2. NORMALIZE (TPM assumed; use vst/rlog if raw counts)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# If already log2(TPM+1), skip normalization:" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "norm_counts <- log2(counts_subset[,-1] + 1)  # log2(TPM+1)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "rownames(norm_counts) <- counts_subset\$transcript" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# IF RAW COUNTS, uncomment:" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_subset[,-1]), colData = coldata, design = ~ condition)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# norm_counts <- vst(dds, blind=TRUE)  # Variance stabilizing, memory efficient" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# 3. SELECT TOP VARIABLE GENES (crucial for 10GB+ large sample sets)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do topvar=`grep "${plot}" ${input} | cut -f5 | sort | uniq`; echo "top_var <- ${topvar}  # Adjust: 1000-5000 usually enough for PCA/heatmap" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "row_vars <- rowVars(as.matrix(norm_counts))" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "top_genes <- order(row_vars, decreasing = TRUE)[1:top_var]" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "norm_top <- norm_counts[top_genes, ]  # e.g. ~2000 x N_samples matrix (manageable)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "cat(\"Top variable genes selected:\", nrow(norm_top), \"\n\")" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# 4. PCA ANALYSIS" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# =============================================================================" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "pca_data <- prcomp(t(norm_top), scale. = TRUE)  # Transpose: samples x genes" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "percentVar <- round(100 * summary(pca_data)\$importance[2, 1:2])" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# PCA plot" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "pca_df <- data.frame(" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  PC1 = pca_data\$x[,1]," >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  PC2 = pca_data\$x[,2]," >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do condition=`grep "${plot}" ${input} | cut -f4 | sort | uniq`; echo "  condition = coldata\$${condition}," >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  name = rownames(coldata)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo ")" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Plot system" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do shapes=`grep "${plot}" ${input} | cut -f6 | sort | uniq`; echo "shape_vals <- c(${shapes})" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = condition)) +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  geom_point(size = 3, alpha = 0.8, stroke = 0.5) +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  xlab(paste0(\"PC1: \", percentVar[1], \"% variance\")) +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  ylab(paste0(\"PC2: \", percentVar[2], \"% variance\")) +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  scale_shape_manual(values = shape_vals) +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  coord_fixed() +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  theme(aspect.ratio = 1) +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  scale_color_viridis_d() +" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  guides(" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "    color = guide_legend(override.aes = list(shape = shape_vals, size = 4))," >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "    shape = guide_legend(override.aes = list(size = 4))" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "  )" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# Save plot" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "ggsave(file.path(dir, \"${outpath_GRADE2081_R}/${plot}\"), plot = pca, width = 8, height = 8, dpi = 300)" >> ${pbs_stem}_${plot}_${thislogdate}.r; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "# File was succesfully written to ${pbs_stem}_${plot}_${thislogdate}.r" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#  Run step" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "#................................................" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo 'echo "## Run R script at" ; date ; echo' >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read plot; do echo "Rscript ${pbs_stem}_${plot}_${thislogdate}.r" >> ${pbs_stem}_${plot}_${thislogdate}.pbs ; done


#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n17 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n15 ${pbs} | head -n1)" ; echo "#                       $(tail -n14 ${pbs} | head -n1)" ; echo "#                       $(tail -n13 ${pbs} | head -n1)" ; echo "#                       $(tail -n11 ${pbs} | head -n1)" ; echo "#                       $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n9 ${pbs} | head -n1)" ; echo "#                       $(tail -n8 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n5 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n3 ${pbs} | head -n1)" ; echo "#                       $(tail -n2 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs} | head -n1)" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including GRADE2 step 081 jobs) at
qstat -u "$user"

# This will remove $VARNAMES from output file with the actual $VARVALUE
# allowing for easily retracing commands
sed -i 's,${input},'"${input}"',g' "$logfile"
sed -i 's,${pbs_stem},'"${pbs_stem}"',g' "$logfile"
sed -i 's,${email},'"${email}"',g' "$logfile"
sed -i 's,${mem},'"${mem}"',g' "$logfile"
sed -i 's,${ncpus},'"${ncpus}"',g' "$logfile"
sed -i 's,${walltime},'"${walltime}"',g' "$logfile"
sed -i 's,${human_thislogdate},'"${human_thislogdate}"',g' "$logfile"
sed -i 's,${thislogdate},'"${thislogdate}"',g' "$logfile"
sed -i 's,${user},'"${user}"',g' "$logfile"
sed -i 's,${module_R},'"${module_R}"',g' "$logfile"
sed -i 's,${outpath_GRADE2081_R},'"${outpath_GRADE2081_R}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v