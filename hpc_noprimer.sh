set -e

export DIR_AMIN=/hpc/cog_bioinf/ridder/users/aallahyar/

# Run specific, should not be hard coded (TODO)
export EXPNAME=$1 	#'LVR-HS2-96x'
export FILE_FASTQ=$2 	#'/home/cog/rstraver/Workspace/pymc4c_hpc/data/NPS_Files/NPS_LVR-HS2-96x_--000010--_Best.fastq'
export FILE_REF=$3 	#$DIR_AMIN/Dataset/Genome_Assembly/Mus_Musculus/mm9/chrAll.fa

# Where the tool/scripts are located
export DIR_TOOL="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export MC4CTOOL=$DIR_TOOL/mc4c.py

# External tools used in the pipeline
export BOWTIE=$DIR_AMIN/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2
export BWA=$DIR_AMIN/My_Works/Useful_Sample_Codes/BWA/bwa/bwa

# Directories for input/output
export DIR_DATA=$DIR_TOOL/data
export DIR_LOG=$DIR_DATA/log
export DIR_NODES=$DIR_TOOL/hpc

# File naming for job input.ouput
export FILE_OUT=$DIR_DATA/$EXPNAME
export FILE_PRIMERFA=$DIR_TOOL/data/primer_${EXPNAME}.fa
export FILE_DATAINFO=$DIR_DATA/ini63/${EXPNAME}.tsv

# Setup for submission variables
QSUBVARS="-V -e $DIR_LOG -o $DIR_LOG"
HOLD_ID=-1

# Load necessary modules for analyses
module load python

# Local: Make the primer fa
#python $MC4CTOOL makeprimerfa $FILE_DATAINFO $FILE_PRIMERFA

# Local: Determine amount of reads in the original fastq file
FASTQWCL=`wc -l $FILE_FASTQ | tail -n 1 |  awk '{print $1}'`
READSPERFILE=16384
export LINESPERFILE=$(($READSPERFILE*4))
export SPLITFILESNUM=$((($FASTQWCL+$LINESPERFILE-1)/$LINESPERFILE))

HOLD_ID_LIST=">"

#echo $(($FASTQWCL/4)) Reads found in $FILE_FASTQ
#echo Becomes $SPLITFILESNUM files

# Single: Split fastq data over several even sized files,
# 	rename reads to unique numbers
HOLD_ID=`qsub $QSUBVARS $DIR_NODES/splitfq.sh | awk '{print $3}'`
HOLD_ID_LIST=$HOLD_ID_LIST"splitfq:$HOLD_ID;"

# Array: Map primers to reads
#HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID -t 1:$SPLITFILESNUM $DIR_NODES/bowtie.sh | awk '{print $3}' | cut -f1 -d.`
#HOLD_ID_LIST=$HOLD_ID_LIST"bowtie:$HOLD_ID;"

# Array: Split reads by primers
#HOLD_ID=`qsub $QSUBVARS -hold_jid_ad $HOLD_ID -t 1:$SPLITFILESNUM $DIR_NODES/splitpr.sh | awk '{print $3}' | cut -f1 -d.`
#HOLD_ID_LIST=$HOLD_ID_LIST"splitpr:$HOLD_ID;"

# Array: Split reads by restriction sites
HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID -t 1:$SPLITFILESNUM $DIR_NODES/splitre_noprimer.sh | awk '{print $3}' | cut -f1 -d.`
HOLD_ID_LIST=$HOLD_ID_LIST"splitre:$HOLD_ID;"

# Single: Merge split reads to one zipped file
HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID $DIR_NODES/mergere.sh | awk '{print $3}'`
HOLD_ID_LIST=$HOLD_ID_LIST"mergere:$HOLD_ID;"

# Single: Map split reads to reference genome
HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID $DIR_NODES/bwa.sh | awk '{print $3}'`
HOLD_ID_LIST=$HOLD_ID_LIST"bwa:$HOLD_ID;"

echo Jobs submitted: $HOLD_ID_LIST
