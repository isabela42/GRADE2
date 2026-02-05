#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Mainá Bitar's 'GRADE (Basic Rnaseq Analysis IN) PBS'
  - Isabela Almeida's 'HyDRA (Hybrid de novo RNA assembly) pipeline'
Created on Jun 18, 2024
Last modified on Feb 05, 2026
Version: ${version}

Description: Write and submit PBS jobs for step 055 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_step055-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources baseline: -m 1 -c 1 -w "20:00:00" #max recources usage with over 7k samples to compile

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            path/from/working/dir/to/grade054_quant_RSEM_DATE/rquant-stem_sampleID.isoforms.results

                            Col2:
                            sampleID
                            as it appears in Col1

                            Col3:
                            expNAME
                            e.g. ovarianRNAseq

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
        ?) echo script usage: bash ipda_GRADE2_step055-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_GRADE2_step055-to-pbs_${thislogdate}.txt

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
outpath_GRADE2055_Bash="grade055_rquant_Bash_${thislogdate}"

## Create output directories
mkdir -p ${outpath_GRADE2055_Bash}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_GRADE2_step055-to-pbs.sh"
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
echo "## Executing bash ipda_GRADE2_step055-to-pbs.sh"
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
cut -f3 ${input} | sort | uniq | while read exp; do echo "#!/bin/sh" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "##########################################################################" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Script:  ${pbs_stem}_${exp}_${thislogdate}.pbs" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Version: v01" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Email:   ${email}" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "##########################################################################" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done

## Write PBS directives
cut -f3 ${input} | sort | uniq | while read exp; do echo "#PBS -N ${pbs_stem}_${exp}_${thislogdate}" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#PBS -r n" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#PBS -m abe" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#PBS -M ${email}" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done

## Write directory setting
cut -f3 ${input} | sort | uniq | while read exp; do echo "#................................................" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Set main working directory" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#................................................" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "## Change to main directory" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done

## Write PBS command lines
cut -f3 ${input} | sort | uniq | while read exp; do echo "#................................................" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#  Run step" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "#................................................" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo 'echo "## Get Expected Counts tables at" ; date ; echo' >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "grep \"${exp}\" ${input} | cut -f2 | tr '\n' '\t' | sed 's/^/transcript\tgene\t/' > ${outpath_GRADE2055_Bash}/${exp}-out-expcnt.tmp" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "sort -k1 \$(grep \"${exp}\" ${input} | cut -f1 | head -n1) | cut -f1-2 | grep -v \"transcript_id\" > ${outpath_GRADE2055_Bash}/${exp}-transcripts-expcnt.tmp" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "grep \"${exp}\" ${input} | while IFS=\$'\t' read -r filepath sampleid exp; do [[ -z \"\${filepath}\" || ! -f \"\${filepath}\" ]] && continue ; sort -k1 \${filepath} | grep -v \"transcript_id\" | cut -f5 > ${outpath_GRADE2055_Bash}/${exp}-tmp-expcnt_\${sampleid} ; done" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "n=100 ; current="${outpath_GRADE2055_Bash}/${exp}-transcripts-expcnt.tmp" ; files=(\$(ls ${outpath_GRADE2055_Bash}/${exp}-tmp-expcnt_* 2>/dev/null )) ; for ((i=0; i<\${#files[@]}; i+=\$n)); do paste -d\$'\t' \$current "\${files[@]:i:\$n}" > ${outpath_GRADE2055_Bash}/${exp}-data-expcnt.tmp.\$\$ ; mv ${outpath_GRADE2055_Bash}/${exp}-data-expcnt.tmp.\$\$ \$current ; done" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "cat ${outpath_GRADE2055_Bash}/${exp}-out-expcnt.tmp ${outpath_GRADE2055_Bash}/${exp}-transcripts-expcnt.tmp > ${outpath_GRADE2055_Bash}/RSEM-${exp}_expcnt.tsv && rm -f ${outpath_GRADE2055_Bash}/${exp}-out-expcnt.tmp ${outpath_GRADE2055_Bash}/${exp}-transcripts-expcnt.tmp ${outpath_GRADE2055_Bash}/${exp}-tmp-expcnt_* " >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done

cut -f3 ${input} | sort | uniq | while read exp; do echo 'echo "## Get TPM tables at" ; date ; echo' >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "grep \"${exp}\" ${input} | cut -f2 | tr '\n' '\t' | sed 's/^/transcript\tgene\t/' > ${outpath_GRADE2055_Bash}/${exp}-out-tpm.tmp" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "sort -k1 \$(grep \"${exp}\" ${input} | cut -f1 | head -n1) | cut -f1-2 | grep -v \"transcript_id\" > ${outpath_GRADE2055_Bash}/${exp}-transcripts-tpm.tmp" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "grep \"${exp}\" ${input} | while IFS=\$'\t' read -r filepath sampleid exp; do [[ -z \"\${filepath}\" || ! -f \"\${filepath}\" ]] && continue ; sort -k1 \${filepath} | grep -v \"transcript_id\" | cut -f6 > ${outpath_GRADE2055_Bash}/${exp}-tmp-tpm_\${sampleid} ; done" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "n=100 ; current="${outpath_GRADE2055_Bash}/${exp}-transcripts-tpm.tmp" ; files=(\$(ls ${outpath_GRADE2055_Bash}/${exp}-tmp-tpm_* 2>/dev/null )) ; for ((i=0; i<\${#files[@]}; i+=\$n)); do paste -d\$'\t' \$current "\${files[@]:i:\$n}" > ${outpath_GRADE2055_Bash}/${exp}-data-tpm.tmp.\$\$ ; mv ${outpath_GRADE2055_Bash}/${exp}-data-tpm.tmp.\$\$ \$current ; done" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "cat ${outpath_GRADE2055_Bash}/${exp}-out-tpm.tmp ${outpath_GRADE2055_Bash}/${exp}-transcripts-tpm.tmp > ${outpath_GRADE2055_Bash}/RSEM-${exp}_tpm.tsv && rm -f ${outpath_GRADE2055_Bash}/${exp}-out-tpm.tmp ${outpath_GRADE2055_Bash}/${exp}-transcripts-tpm.tmp ${outpath_GRADE2055_Bash}/${exp}-tmp-tpm_* " >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done

cut -f3 ${input} | sort | uniq | while read exp; do echo 'echo "## Get FPKM tables at" ; date ; echo' >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "grep \"${exp}\" ${input} | cut -f2 | tr '\n' '\t' | sed 's/^/transcript\tgene\t/' > ${outpath_GRADE2055_Bash}/${exp}-out-fpkm.tmp" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "sort -k1 \$(grep \"${exp}\" ${input} | cut -f1 | head -n1) | cut -f1-2 | grep -v \"transcript_id\" > ${outpath_GRADE2055_Bash}/${exp}-transcripts-fpkm.tmp" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "grep \"${exp}\" ${input} | while IFS=\$'\t' read -r filepath sampleid exp; do [[ -z \"\${filepath}\" || ! -f \"\${filepath}\" ]] && continue ; sort -k1 \${filepath} | grep -v \"transcript_id\" | cut -f7 > ${outpath_GRADE2055_Bash}/${exp}-tmp-fpkm_\${sampleid} ; done" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "n=100 ; current="${outpath_GRADE2055_Bash}/${exp}-transcripts-fpkm.tmp" ; files=(\$(ls ${outpath_GRADE2055_Bash}/${exp}-tmp-fpkm_* 2>/dev/null )) ; for ((i=0; i<\${#files[@]}; i+=\$n)); do paste -d\$'\t' \$current "\${files[@]:i:\$n}" > ${outpath_GRADE2055_Bash}/${exp}-data-fpkm.tmp.\$\$ ; mv ${outpath_GRADE2055_Bash}/${exp}-data-fpkm.tmp.\$\$ \$current ; done" >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done
cut -f3 ${input} | sort | uniq | while read exp; do echo "cat ${outpath_GRADE2055_Bash}/${exp}-out-fpkm.tmp ${outpath_GRADE2055_Bash}/${exp}-transcripts-fpkm.tmp > ${outpath_GRADE2055_Bash}/RSEM-${exp}_fpkm.tsv && rm -f ${outpath_GRADE2055_Bash}/${exp}-out-fpkm.tmp ${outpath_GRADE2055_Bash}/${exp}-transcripts-fpkm.tmp ${outpath_GRADE2055_Bash}/${exp}-tmp-fpkm_* " >> ${pbs_stem}_${exp}_${thislogdate}.pbs ; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n17 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n15 ${pbs} | head -n1)" ; echo "#                       $(tail -n14 ${pbs} | head -n1)" ; echo "#                       $(tail -n13 ${pbs} | head -n1)" ; echo "#                       $(tail -n11 ${pbs} | head -n1)" ; echo "#                       $(tail -n10 ${pbs} | head -n1)" ; echo "#                       $(tail -n9 ${pbs} | head -n1)" ; echo "#                       $(tail -n8 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n5 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#                       $(tail -n3 ${pbs} | head -n1)" ; echo "#                       $(tail -n2 ${pbs} | head -n1)" ; echo "#                       $(tail -n1 ${pbs} | head -n1)" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including GRADE2 step 055 jobs) at
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
sed -i 's,${outpath_GRADE2055_Bash},'"${outpath_GRADE2055_Bash}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v