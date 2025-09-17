#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Mainá Bitar's 'GRADE (Basic Rnaseq Analysis IN) PBS'
  - Isabela Almeida's 'HyDRA (Hybrid de novo RNA assembly) pipeline'
Created on Jun 18, 2024
Last modified on September 16, 2025
Version: ${version}

Description: Write and submit PBS jobs for step 091 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_step091-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources baseline: -m 2 -c 1 -w "10:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            path/from/working/dir/to/grade081_quant_RSEM_DATE/

-p <PBS stem>               Stem for PBS file names
-e <email>                  Email for PBS job
-m <INT>                    Memory INT required for PBS job in GB
-c <INT>                    Number of CPUS required for PBS job
-w <HH:MM:SS>               Clock walltime required for PBS job

## Output:

logfile.txt                 Logfile with commands executed and date
PBS files                   PBS files created

Pipeline description:

#   000 Index building (1Kallisto, 2RSEM, 3STAR)
#   010 Quality check raw files (1FastQC, 2MultiQC)
#   020 Trim reads of adapters (1Trimmomatic)
#   030 Quality check raw files (1FastQC, 2MultiQC)
#   040 Quantify reads (1Kallisto)
#   050 Create Kallisto count tables (1Kallisto)
#   060 Alignment (1STAR)
#   070 Process alignment (1SAMtools, 2NovoSort)
#   080 Quantify reads (1RSEM)
#-->090 Create RSEM count tables (1RSEM)
#   100 Differential Expression Analysis (1EdgeR)

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
        ?) echo script usage: bash ipda_GRADE2_step091-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_GRADE2_step091-to-pbs_${thislogdate}.txt

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_GRADE2_step091-to-pbs.sh"
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
echo "## Executing bash ipda_GRADE2_step091-to-pbs.sh"
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
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#!/bin/sh" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "##########################################################################" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Script:  ${pbs_stem}_${folder_date}_${thislogdate}.pbs" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Version: v01" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Email:   ${email}" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "##########################################################################" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#PBS -N ${pbs_stem}_${thislogdate}" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#PBS -r n" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#PBS -m abe" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#PBS -M ${email}" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#................................................" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Set main working directory" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#................................................" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "## Change to main directory" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#................................................" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#  Run step" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "#................................................" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo "## Get expected counts tables at" ; date ; echo' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "ls -d ${path_input}/rquant*isoforms.results | cut -d"." -f1 | while read folder; do seed=\`echo \${folder} | rev | cut -d\"/\" -f1 | rev | cut -d\"_\" -f4-\` ; echo -e "Transcript\t${seed}" >> \${folder}_ExpCounts; cut -f1,5 \${folder}.isoforms.results | tail -n+2 >> \${folder}_ExpCounts; done" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "first=\`ls -d ${path_input}/rquant*ExpCounts | head -n2 | head -n1\`; second=\`ls -d ${path_input}/rquant*ExpCounts | head -n2 | tail -n1\`; n=\`ls -d ${path_input}/rquant*ExpCounts | wc -l\`; ns=\`echo \"\${n}-2\" | bc\`; np=\`echo \"\${n}+1\" | bc\`" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "join -t\$'\t' \${first} \${second} >> ${path_input}/Seed_ExpCounts" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "ls -d ${path_input}/rquant*ExpCounts | tail -n\${ns} | while read folder; do mv ${path_input}/Seed_ExpCounts ${path_input}/Active_ExpCounts; join -t\$'\t' ${path_input}/Active_ExpCounts \${folder} >> ${path_input}/Seed_ExpCounts; done" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "mv ${path_input}/Seed_ExpCounts ${path_input}/Rsem_all_ExpCounts ; rm -f ${path_input}/Active_ExpCounts" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo "## Sneak peak at expected counts table at" ; date ; echo' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "head ${path_input}/Rsem_all_ExpCounts" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo "## Get TPM tables at" ; date ; echo' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "ls -d ${path_input}/rquant*isoforms.results | cut -d"." -f1 | while read folder; do seed=\`echo \${folder} | rev | cut -d\"/\" -f1 | rev | cut -d\"_\" -f4-\` ; echo -e "Transcript\t${seed}" >> \${folder}_TPM; cut -f1,6 \${folder}.isoforms.results | tail -n+2 >> \${folder}_TPM; done" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "first=\`ls -d ${path_input}/rquant*TPM | head -n2 | head -n1\`; second=\`ls -d ${path_input}/rquant*TPM | head -n2 | tail -n1\`; n=\`ls -d ${path_input}/rquant*TPM | wc -l\`; ns=\`echo \"\${n}-2\" | bc\`; np=\`echo \"\${n}+1\" | bc\`" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "join -t\$'\t' \${first} \${second} >> ${path_input}/Seed_TPM" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "ls -d ${path_input}/rquant*TPM | tail -n\${ns} | while read folder; do mv ${path_input}/Seed_TPM ${path_input}/Active_TPM; join -t\$'\t' ${path_input}/Active_TPM \${folder} >> ${path_input}/Seed_TPM; done" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "mv ${path_input}/Seed_TPM ${path_input}/Rsem_all_TPM ; rm -f ${path_input}/Active_TPM" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo "## Sneak peak at TPM counts table at" ; date ; echo' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "head ${path_input}/Rsem_all_TPM" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo "## Get FPKM tables at" ; date ; echo' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "ls -d ${path_input}/rquant*isoforms.results | cut -d"." -f1 | while read folder; do seed=\`echo \${folder} | rev | cut -d\"/\" -f1 | rev | cut -d\"_\" -f4-\` ; echo -e "Transcript\t${seed}" >> \${folder}_FPKM; cut -f1,7 \${folder}.isoforms.results | tail -n+2 >> \${folder}_FPKM; done" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "first=\`ls -d ${path_input}/rquant*FPKM | head -n2 | head -n1\`; second=\`ls -d ${path_input}/rquant*FPKM | head -n2 | tail -n1\`; n=\`ls -d ${path_input}/rquant*FPKM | wc -l\`; ns=\`echo \"\${n}-2\" | bc\`; np=\`echo \"\${n}+1\" | bc\`" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "join -t\$'\t' \${first} \${second} >> ${path_input}/Seed_FPKM" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "ls -d ${path_input}/rquant*FPKM | tail -n\${ns} | while read folder; do mv ${path_input}/Seed_FPKM ${path_input}/Active_FPKM; join -t\$'\t' ${path_input}/Active_FPKM \${folder} >> ${path_input}/Seed_FPKM; done" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "mv ${path_input}/Seed_FPKM ${path_input}/Rsem_all_FPKM ; rm -f ${path_input}/Active_FPKM" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo 'echo "## Sneak peak at FPKM counts table at" ; date ; echo' >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done
cut -f1 ${input} | sort | uniq | while read path_input; do folder_date=`echo ${path_input} | rev | cut -d"_" -f1 | rev`; echo "head ${path_input}/Rsem_all_FPKM" >> ${pbs_stem}_${folder_date}_${thislogdate}.pbs ; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n26 ${pbs} | head -n1)" ; echo "#                       $(tail -n25 ${pbs} | head -n1)" ; echo "#                       $(tail -n28 ${pbs} | head -n1)" ; echo "#                       $(tail -n27 ${pbs} | head -n1)" ; echo "#                       $(tail -n26 ${pbs} | head -n1)" ; echo "#                       $(tail -n25 ${pbs} | head -n1)" ; echo "#                       $(tail -n24 ${pbs} | head -n1)" ; echo "#                       $(tail -n18 ${pbs} | head -n1)" ; echo "#                       $(tail -n17 ${pbs} | head -n1)" ; echo "#                       $(tail -n16 ${pbs} | head -n1)" ; echo "#                       $(tail -n15 ${pbs} | head -n1)" ; echo "#                       $(tail -n14 ${pbs} | head -n1)" ; echo "#                       $(tail -n8 ${pbs} | head -n1)" ; echo "#                       $(tail -n7 ${pbs} | head -n1)" ; echo "#                       $(tail -n6 ${pbs} | head -n1)" ; echo "#                       $(tail -n5 ${pbs} | head -n1)" ; echo "#                       $(tail -n4 ${pbs} | head -n1)" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including GRADE2 step 091 jobs) at
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
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v