#!/bin/bash

version="1.0.0"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Larissa Cassiano's script
Created on Sep 11, 2026
Last modified on Sep 24, 2026
Version: ${version}

Description: Write and submit PBS jobs for step 071 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_step071-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources used for pipeline in-house: -m 1 -c 1 -w "01:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            experiment-stem
                            Note: Can be subgroups present in the same counts files

                            Col2:
                            /path/from/working/dir/to/grade055_rquant_Bash_DATE/RSEM-stem
                            Note: path to stem of TPM and EXPCNT TSV files

                            Col3:
                            /path/from/working/dir/to/annotation.tsv
                            Note: annotation file with GENE_ID|TRANSCRIPT_ID|GENE_SYMBOL|TRANSCRIPT_SYMBOL
                            where IDs match those of inputexp and inputtpm (header false)

                            Col4:
                            groups
                            Note: Comma separated group labels based on order of EXP/TPM counts file, e.g.:
                            \"CTRL,ASO1,ASO2,CTRL,ASO1,ASO2,CTRL,ASO1,ASO2\"

                            Col5:
                            col number INT
                            Note: Column numbers from inputexp and inputtpm to use.
                            Always include 1,2 (transcript,gene) followed by the other column numbers.

                            Col6:
                            contrast groups
                            Note: Pairwise comma separated raw contrast groups for differential expression analysis,
                            e.g.: \"ASO1_vs_CTRL,ASO2_vs_CTRL\"

                            Col7:
                            merge groups
                            Note: Merging rule for group labels, e.g.: \"ASO1,ASO2=ASO\"

                            Col8:
                            merge contrast groups
                            Note: Pairwise comma separated merged contrast groups for differential 
                            expression analysis, e.g.: \"ASO_vs_CTRL\"

                            Col9:
                            genesets
                            Note: Highlight inbuilt gene sets list on plots (use \"NA\" for none,
                            options: HR,FA,BER,NHEJ,Cell_Cycle,top10DEG)

                            Col10:
                            label for plots

                            Col11:
                            /path/from/working/dir/to/msigdb_collections/human
                            Note: Input GSEA MSigDB collections folder. Download from
                            https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

                            Col12:
                            /path/from/working/dir/to/GRADE2/scripts/ipda_grade2_step071.r

                            Col13:
                            /path/from/working/dir/to/GRADE2/scripts/ipda_grade2_rfunctions.r

                            It does not matter if same stem 
                            appears more than once on this input file.

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
#   050 Align (1STAR, 2SAMtools, 3NovoSort) and quantify reads (4RSEM, 5BASH count tables)
#   060 PSeudo align and quantify reads at isoform level (1Salmon, 2BASH count tables)
#-->070 Differential Expression Analysis (1EdgeR)
#   080 Plot counts and draw metrics (1pca+heatmap+boxplot, 2metrics, 3id plots)

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
        ?) echo script usage: bash ipda_grade2_step071-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_GRADE2_step071-to-pbs_${thislogdate}.txt

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
outpath_GRADE2071_plots="grade071_DE-pathway-plots_R_${thislogdate}"

## Create output directories
mkdir -p ${outpath_GRADE2071_plots}

#................................................
#  Required modules, softwares and libraries
#................................................

# R 4.5.0
module_R="R/4.5.0"

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_grade2_step071-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input files:                 ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${outpath_GRADE2071_plots}"
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
echo "## Executing bash ipda_grade2_step071-to-pbs.sh"
echo "## This execution PID: ${pid}"
echo
echo "## Given inputs:"
echo
echo "## Input files:                 ${input}"
echo "## PBS stem:                    ${pbs_stem}"
echo "## Email for PBS notifications: ${email}"
echo "## PBS job memory required:     ${mem}"
echo "## PBS job NCPUS required:      ${ncpus}"
echo "## PBS job walltime required:   ${walltime}"
echo
echo "## Outputs created:"
echo
echo "## Output files saved to:       ${outpath_GRADE2071_plots}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read stem; do echo "#!/bin/sh" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "##########################################################################" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Script:  ${pbs_stem}_${stem}_${thislogdate}.pbs" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Version: v01" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Email:   ${email}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "##########################################################################" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read stem; do echo "#PBS -N ${pbs_stem}_${stem}_${thislogdate}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#PBS -r n" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#PBS -m abe" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#PBS -M ${email}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Set main working directory" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "## Change to main directory" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "module load ${module_R}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done


## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#  Run step" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do echo 'echo "## Run edgeR and plot at" ; date ; echo' >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do inputcnt=`grep "${stem}" ${input} | cut -f2 | sort | uniq`; inputanno=`grep "${stem}" ${input} | cut -f3 | sort | uniq`; groups=`grep "${stem}" ${input} | cut -f4 | sort | uniq`; colnumns=`grep "${stem}" ${input} | cut -f5 | sort | uniq`; contrasts=`grep "${stem}" ${input} | cut -f6 | sort | uniq`; merge=`grep "${stem}" ${input} | cut -f7 | sort | uniq`; mergecontrasts=`grep "${stem}" ${input} | cut -f8 | sort | uniq`; genesets=`grep "${stem}" ${input} | cut -f9 | sort | uniq`; label=`grep "${stem}" ${input} | cut -f10 | sort | uniq`; inputgmt=`grep "${stem}" ${input} | cut -f11 | sort | uniq`; rscript=`grep "${stem}" ${input} | cut -f12 | sort | uniq`; rfunctions=`grep "${stem}" ${input} | cut -f13 | sort | uniq`; echo "Rscript ${rscript} --inputanno ${inputanno} --inputexp ${inputcnt}_isoforms-expcnt.tsv --inputtpm ${inputcnt}_isoforms-tpm.tsv --level \"transcript\" --groups ${groups} --colnumns ${colnumns} --contrasts ${contrasts} --merge ${merge} --mergecontrasts ${mergecontrasts} --genesets ${genesets} --plotlabel ${label} --inputgmt ${inputgmt} --outdir ${outpath_GRADE2071_plots}/ --outstem \"${stem}\" --function ${rfunctions}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read stem; do inputcnt=`grep "${stem}" ${input} | cut -f2 | sort | uniq`; inputanno=`grep "${stem}" ${input} | cut -f3 | sort | uniq`; groups=`grep "${stem}" ${input} | cut -f4 | sort | uniq`; colnumns=`grep "${stem}" ${input} | cut -f5 | sort | uniq`; contrasts=`grep "${stem}" ${input} | cut -f6 | sort | uniq`; merge=`grep "${stem}" ${input} | cut -f7 | sort | uniq`; mergecontrasts=`grep "${stem}" ${input} | cut -f8 | sort | uniq`; genesets=`grep "${stem}" ${input} | cut -f9 | sort | uniq`; label=`grep "${stem}" ${input} | cut -f10 | sort | uniq`; inputgmt=`grep "${stem}" ${input} | cut -f11 | sort | uniq`;rscript=`grep "${stem}" ${input} | cut -f12 | sort | uniq`; rfunctions=`grep "${stem}" ${input} | cut -f13 | sort | uniq`; echo "Rscript ${rscript} --inputanno ${inputanno} --inputexp ${inputcnt}_isoforms-expcnt.tsv --inputtpm ${inputcnt}_genes-tpm.tsv --level \"gene\" --groups ${groups} --colnumns ${colnumns} --contrasts ${contrasts} --merge ${merge} --mergecontrasts ${mergecontrasts} --genesets ${genesets} --plotlabel ${label} --inputgmt ${inputgmt} --outdir ${outpath_GRADE2071_plots}/ --outstem \"${stem}\" --function ${rfunctions}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n2 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs})" ;echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including GRADE2 step 071 jobs) at
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
sed -i 's,${outpath_GRADE2071_plots},'"${outpath_GRADE2071_plots}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v