#!/bin/bash

version="1.0.02"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Mainá Bitar's 'GRADE (Basic Rnaseq Analysis IN) PBS'
  - Isabela Almeida's 'HyDRA (Hybrid de novo RNA assembly) pipeline'
Created on May 27, 2024
Last modified on November 06, 2025
Version: ${version}

Description: Write and submit PBS jobs for step 012 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_step012-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources baseline: -m 1 -c 1 -w "01:00:00"

## Input:

-i <path/to/input/folder>   path/from/working/dir/to/GRADE2_step011_QC-raw-files_FastQC_DATE/
                            This folder should contain the HTML FastQC from step011 files.
-p <PBS stem>               Stem for PBS file names
-e <email>                  Email for PBS job
-m <INT>                    Memory INT required for PBS job in GB
-c <INT>                    Number of CPUS required for PBS job
-w <HH:MM:SS>               Clock walltime required for PBS job

## Output:

logfile.txt                 Logfile with commands executed and date
PBS files                   PBS files created

Pipeline description:

#   000 Index building (0BedTools, 1Kallisto, 2RSEM, 3STAR)
#-->010 Quality check raw files (1FastQC, 2MultiQC)
#   020 Trim reads of adapters (1Trimmomatic)
#   030 Quality check raw files (1FastQC, 2MultiQC)
#   040 Quantify reads (1Kallisto)
#   050 Create Kallisto count tables (1Kallisto)
#   060 Alignment (1STAR)
#   070 Process alignment (1SAMtools, 2NovoSort)
#   080 Quantify reads (1RSEM)
#   090 Create RSEM/Kallisto count tables (1Kallisto-RSEM)
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
        ?) echo script usage: bash ipda_GRADE2_step012-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_GRADE2_step012-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# MultiQC:
# <https://multiqc.info/docs/getting_started/installation/>
multiqc="/working/lab_julietF/isabelaA/tools/anaconda3/bin/multiqc"

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
outpath_GRADE2012_MultiQC="grade012_qc-raw_MultiQC_${thislogdate}"

## Create output directories
mkdir -p ${outpath_GRADE2012_MultiQC}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_GRADE2_step012-to-pbs.sh"
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
echo "## Output files saved to:       ${outpath_GRADE2012_MultiQC}"
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
echo "## Executing bash ipda_GRADE2_step012-to-pbs.sh"
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
echo "## Output files saved to:       ${outpath_GRADE2012_MultiQC}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
echo "#!/bin/sh" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo "##########################################################################" >> ${pbs_stem}_${thislogdate}.pbs
echo "#" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Script:  ${pbs_stem}_${thislogdate}.pbs" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Version: v01" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Email:   ${email}" >> ${pbs_stem}_${thislogdate}.pbs
echo "#" >> ${pbs_stem}_${thislogdate}.pbs
echo "##########################################################################" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write PBS directives
echo "#PBS -N ${pbs_stem}_${thislogdate}" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -r n" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -m abe" >> ${pbs_stem}_${thislogdate}.pbs
echo "#PBS -M ${email}" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write directory setting
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Set main working directory" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo "## Change to main directory" >> ${pbs_stem}_${thislogdate}.pbs
echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write load modules
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo "echo \"${multiqc}\"" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs

## Write PBS command lines
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "#  Run step" >> ${pbs_stem}_${thislogdate}.pbs
echo "#................................................" >> ${pbs_stem}_${thislogdate}.pbs
echo "" >> ${pbs_stem}_${thislogdate}.pbs
echo 'echo "## Run MultiQC at" ; date ; echo' >> ${pbs_stem}_${thislogdate}.pbs
echo "${multiqc} ${input} -o ${outpath_GRADE2012_MultiQC}" >> ${pbs_stem}_${thislogdate}.pbs

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (includingGRADE2 step 012 jobs) at
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
sed -i 's,${multiqc},'"${multiqc}"',g' "$logfile"
sed -i 's,${outpath_GRADE2012_MultiQC},'"${outpath_GRADE2012_MultiQC}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v