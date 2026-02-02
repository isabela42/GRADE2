#!/bin/bash

version="1.0.02"
usage(){
echo "
Written by Isabela Almeida
Based on
  - Mainá Bitar's 'GRADE (Basic Rnaseq Analysis IN) PBS'
  - Isabela Almeida's 'HyDRA (Hybrid de novo RNA assembly) pipeline'
Created on Jun 04, 2024
Last modified on Feb 02, 2026
Version: ${version}

Description: Write and submit PBS jobs for step 051 of the
GRADE2 PBS 2.0 pipeline (General RNAseq Analysis for Differential Expression version 2).

Usage: bash ipda_GRADE2_step051-to-pbs.sh -i "path/to/input/files" -p "PBS stem" -e "email" -m INT -c INT -w "HH:MM:SS"

Resources baseline: -m 40 -c 6 -w "10:00:00"

## Input:

-i <path/to/input/files>    Input TSV files including path FROM working
                            directory. This TSV file should contain:
                            
                            Col1:
                            path/from/working/dir/to/grade021_trim-adapters_Trimmomatic_DATE/adapter-trimmed_stem
                            of both _R1.f* and _R2.f* files in individual lines
                            and no full stops.

                            Col2:
                            path/to/star-index

                            Col3:
                            path/from/working/dir/to/file.gtf
                            of user desired transcriptome (e.g. gencode +
                            custom transcriptome in a single gtf file
                            Note: must be the same used to create the index

                            Col4:
                            single|paired
                            Flag for single-end or paired-end reads. If Col4 is single, it will run with F1 only.

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
        ?) echo script usage: bash ipda_GRADE2_step051-to-pbs.sh -i path/to/input/files -p PBS stem -e email -m INT -c INT -w "HH:MM:SS" >&2
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
logfile=logfile_ipda_GRADE2_step051-to-pbs_${thislogdate}.txt

#................................................
#  Required modules, softwares and libraries
#................................................

# STAR 2.7.1a:
module_star=STAR/2.7.1a

#................................................
#  Set and create output path
#................................................

## Set stem for output directories
outpath_GRADE2051_STAR="grade051_alignment_STAR_${thislogdate}"

## Create output directories
mkdir -p ${outpath_GRADE2051_STAR}

#................................................
#  Print Execution info to user
#................................................

date
echo "## Executing bash ipda_GRADE2_step051-to-pbs.sh"
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
echo "## Output files saved to:       ${outpath_GRADE2051_STAR}"
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
echo "## Executing bash ipda_GRADE2_step051-to-pbs.sh"
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
echo "## Output files saved to:       ${outpath_GRADE2051_STAR}"
echo "## This is logfile:             ${logfile}"

set -v

#................................................
#  Create PBS files
#................................................

## Write PBS header
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#!/bin/bash" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Script:  ${pbs_stem}_${file}_${thislogdate}.pbs" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Author:  Isabela Almeida" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Created: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Updated: ${human_thislogdate} at QIMR Berghofer (Brisbane, Australia) - VSC" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Version: v01" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Email:   ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "##########################################################################" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS directives
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#PBS -N ${pbs_stem}_${file}_${thislogdate}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#PBS -r n" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#PBS -l mem=${mem}GB,walltime=${walltime},ncpus=${ncpus}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#PBS -m ae" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#PBS -M ${email}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write directory setting
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Set main working directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "## Change to main directory" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo 'cd ${PBS_O_WORKDIR}' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo 'echo ; echo "WARNING: The main directory for this run was set to ${PBS_O_WORKDIR}"; echo ' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write load modules
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Load Softwares, Libraries and Modules" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "module load ${module_star}" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

## Write PBS command lines
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#  Run step" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "#................................................" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo "" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; echo 'echo "## Run STAR at" ; date ; echo' >> ${pbs_stem}_${file}_${thislogdate}.pbs; done
cut -f1 ${input} | sort | uniq | while read path_file; do file=`echo ${path_file} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2-` ; f1="${path_file}_R1.fq.gz"; f2="${path_file}_R2.fq.gz"; index=`grep -P "${path_file}\t" ${input} | cut -f2`; indexstem=`echo ${index} | rev | cut -d"/" -f1 | rev`; gtf=`grep -P "${path_file}\t" ${input} | cut -f3`; mkdir -p ${outpath_GRADE2051_STAR}/salign-${indexstem}_${file}; flag_reads=`grep -P "${path_file}\t" ${input} | cut -f4`; if [[ $flag_reads = "single" ]]; then reads="${f1}"; else reads="${f1} ${f2}"; fi; echo "STAR --runMode alignReads --readFilesIn ${reads} --readFilesCommand zcat --outFileNamePrefix ${outpath_GRADE2051_STAR}/salign-${indexstem}_${file}/salign-${indexstem}_${file} --genomeDir ${index} --sjdbGTFfile ${gtf} --outSJfilterReads Unique --sjdbOverhang 100 --runThreadN ${ncpus} --genomeLoad NoSharedMemory --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --chimSegmentMin 20 --outSAMattributes All --outSAMstrandField intronMotif --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate" >> ${pbs_stem}_${file}_${thislogdate}.pbs; done

# Star command explained
# genomeDir: path to the directory where genome files are stored (if runMode!=generateGenome) or will be generated (if runMode==generateGenome)

#		--readFilesCommand zcat \			string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
#		--sjdbGTFfile ${gtfDir} \			path to the GTF file with annotations
#		--outSJfilterReads Unique \			string: which reads to consider for collapsed splice junctions output (uniquely mapping reads only)
#		--sjdbOverhang 114 \				int>0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)
#		--twopassMode Basic \				string: 2-pass mapping mode (basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly)
#		--outFilterType BySJout \			string: type of filtering (BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab)
#		--outFilterMultimapNmax 100 \			int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out .
#		--outFilterMismatchNmax 33 \			int: alignment will be output only if it has no more mismatches than this value.
#		--outFilterMatchNminOverLread 0 \		float: same as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads). // int: alignment will be output only if the number of matched bases is higher than or equal to this value.
#		--outFilterMismatchNoverLmax 0.3 \		float: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value.
#		--outFilterScoreMinOverLread 0.3 \		float: same as outFilterScoreMin, but  normalized to read length (sum of mates' lengths for paired-end reads) // int: alignment will be output only if its score is higher than or equal to this value.
#		--limitOutSJcollapsed 1000000 \			int>0: max number of collapsed junctions
#		--limitSjdbInsertNsj 1000000 \			int>=0: maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run
#		--alignSJoverhangMin 8 \			int>0: minimum overhang (i.e. block size) for spliced alignments
#		--alignEndsType EndToEnd \			string: type of read ends alignment (force end-to-end read alignment, do not soft-clip)
#		--alignSJDBoverhangMin 3  \			int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments
#		--alignIntronMin 20 \				maximum intron size (if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins)
#		--winAnchorMultimapNmax 50 \			int>0: max number of loci anchors are allowed to map to
#		--seedSearchStartLmax 12 \			int>0: defines the search start point through the read - the read is split into pieces no longer than this value
#		--chimSegmentMin 20 \				int>=0: minimum length of chimeric segment length (if ==0, no chimeric output)
#		--outSAMattributes All \			string: a string of desired SAM attributes, in the order desired for the output SAM (All = NH HI AS nM NM MD jM jI)
#		--outSAMstrandField intronMotif \		string: Cufflinks-like strand field flag (intronMotif = strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.)
#		--quantMode TranscriptomeSAM \			string(s): types of quantification requested (TranscriptomeSAM = output SAM/BAM alignments to transcriptome into a separate file)
#		--outSAMattrIHstart 0 \				int>=0: start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.
#		--outSAMunmapped Within \			string(s): output of unmapped reads in the SAM format (Within = output unmapped reads within the main SAM file (i.e. Aligned.out.sam))


#................................................
#  Submit PBS jobs
#................................................

## Submit PBS jobs 
ls ${pbs_stem}_*${thislogdate}.pbs | while read pbs; do echo ; echo "#................................................" ; echo "# This is PBS: ${pbs}" ;  echo "#" ; echo "# main command line(s): $(tail -n1 ${pbs})" ; echo "#" ; echo "# now submitting PBS" ; echo "qsub ${pbs}" ; qsub ${pbs} ; echo "#................................................" ; done

date ## Status of all user jobs (including GRADE2 step 051 jobs) at
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
sed -i 's,${module_star},'"${module_star}"',g' "$logfile"
sed -i 's,${outpath_GRADE2051_STAR},'"${outpath_GRADE2051_STAR}"',g' "$logfile"
sed -i 's,${logfile},'"${logfile}"',g' "$logfile"
sed -n -e :a -e '1,3!{P;N;D;};N;ba' $logfile > tmp ; mv tmp $logfile
set +v