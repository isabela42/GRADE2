#!/bin/bash

version="1.0.00"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Mainá Bitar's 'GRADE (Basic Rnaseq Analysis IN) PBS'
  - Isabela Almeida's 'HyDRA (Hybrid de novo RNA assembly) pipeline'
Created on Mar 09, 2024
Last modified on Mar 13, 2026
Version: ${version}

Description: Write and submit PBS jobs for step 082 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_step082-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources baseline: -m 1 -c 1 -w "01:00:00" # scales up depending on # samples and when calculating median

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            /path/from/working/dir/to/grade0*_*quant_Bash_DATE/*-hydra-rnaseq_*.tsv
                            counts file from 042, 082 or 062

                            Col2:
                            /path/from/working/dir/to/target-metadata.tsv
                            Note: this file must contain all target samples ids in col1,
                            i.e. ids from samples that want to be above expression threshold

                            Col3:
                            INT
                            threshold value, e.g. 0.5

                            Col4:
                            INT
                            proportion of targets that must be above col3 threshold

                            Col5:
                            INT
                            proportion of non-targets that must be below col3 threshold
                            
                            Col6:
                            stem-of-filter
                            Note: stem-of-filter can only appear once per input file

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
#   070 Differential Expression Analysis (1EdgeR)
#-->080 Plot counts and draw metrics (1pca, 2metrics)

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
        ?) echo script usage: bash ipda_GRADE2_step082-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_GRADE2_step082-to-pbs_${thislogdate}.txt

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
outpath_GRADE2082_Bash="grade082_metrics_Bash_${thislogdate}"

## Create output directories
mkdir -p ${outpath_GRADE2082_Bash}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_GRADE2_step082-to-pbs.sh"
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
echo "## Executing bash ipda_GRADE2_step082-to-pbs.sh"
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
cut -f6 ${input} | sort | uniq | while read stem; do echo "#!/bin/sh" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "##########################################################################" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Script:  ${pbs_stem}_${stem}_${thislogdate}.pbs" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Version: v01" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Email:   ${email}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "##########################################################################" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done

## Write PBS directives
cut -f6 ${input} | sort | uniq | while read stem; do echo "#PBS -N ${pbs_stem}_${stem}_${thislogdate}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#PBS -r n" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#PBS -m abe" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#PBS -M ${email}" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done

## Write directory setting
cut -f6 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Set main working directory" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "## Change to main directory" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done

## Write PBS command lines
cut -f6 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#  Run step" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "#................................................" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo "" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do echo 'echo "## Filter ids at" ; date ; echo' >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; done
cut -f6 ${input} | sort | uniq | while read stem; do counts=`grep "${stem}" ${input} | cut -f1 | sort | uniq`; metadata=`grep "${stem}" ${input} | cut -f2 | sort | uniq`; thr=`grep "${stem}" ${input} | cut -f3 | sort | uniq`; propTarg=`grep "${stem}" ${input} | cut -f4 | sort | uniq`;  propNTarg=`grep "${stem}" ${input} | cut -f5 | sort | uniq`; if [ "${propTarg}" == "median" ]; then echo "awk -v thr=${thr} -v pn=${propNTarg} 'NR==FNR{targets[\$1]=1;next} FNR==1{for(i=2;i<=NF;i++){if(\$i in targets){tcol[i]=1;nt++}else{ocol[i]=1;no++}};next} {k=0;hit_o=0;for(i in tcol){k++;v[k]=\$i} for(i=1;i<=k;i++)for(j=i+1;j<=k;j++)if(v[i]>v[j]){tmp=v[i];v[i]=v[j];v[j]=tmp} med=(k%2?v[(k+1)/2]:(v[k/2]+v[k/2+1])/2); for(i in ocol) if(\$i>=thr) hit_o++; if(med>=thr && hit_o<=pn*no) print \$1}' ${metadata} ${counts} >> ${outpath_GRADE2082_Bash}/${stem}.txt" >> ${pbs_stem}_${stem}_${thislogdate}.pbs; elif [ "${propTarg}" == "one" ]; then echo "awk -v thr=${thr} -v pt=${propTarg} -v pn=${propNTarg} 'NR==FNR{targets[\$1]=1;next} FNR==1{for(i=2;i<=NF;i++){if(\$i in targets){tcol[i]=1;nt++}else{ocol[i]=1;no++}};next} {hit_t=0;hit_o=0;for(i in tcol) if(\$i>=thr) hit_t++; for(i in ocol) if(\$i>=thr) hit_o++; if(hit_t>=1 && hit_o<=pn*no) print \$1}' ${metadata} ${counts} >> ${outpath_GRADE2082_Bash}/${stem}.txt" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; else echo "awk -v thr=${thr} -v pt=${propTarg} -v pn=${propNTarg} 'NR==FNR{targets[\$1]=1;next} FNR==1{for(i=2;i<=NF;i++){if(\$i in targets){tcol[i]=1;nt++}else{ocol[i]=1;no++}};next} {hit_t=0;hit_o=0;for(i in tcol) if(\$i>=thr) hit_t++; for(i in ocol) if(\$i>=thr) hit_o++; if(hit_t>=pt*nt && hit_o<=pn*no) print \$1}' ${metadata} ${counts} >> ${outpath_GRADE2082_Bash}/${stem}.txt" >> ${pbs_stem}_${stem}_${thislogdate}.pbs ; fi; done

#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n1 ${pbs})"; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including GRADE2 step 082 jobs) at
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
sed -i 's,${outpath_GRADE2082_Bash},'"${outpath_GRADE2082_Bash}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v