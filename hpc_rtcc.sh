
set -e

## Debug mode?
export DEBUG_MODE="set -o xtrace"
$DEBUG_MODE

# Run specific
export FILE_INI=`realpath $1`
export DIR_WORKSPACE=`realpath $2`

# External tools used in the pipeline
export DIR_AMIN=/hpc/cog_bioinf/ridder/users/aallahyar/
export BOWTIE=$DIR_AMIN/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2
export BWA=$DIR_AMIN/My_Works/Useful_Sample_Codes/BWA/bwa/bwa

# Where the tool/scripts are located
export DIR_MC4C="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export MC4CTOOL=$DIR_MC4C/mc4c.py

# Directories for input/output

SOURCE_FASTQ=`awk '/^src_fastq/{$1="";print $0}' $FILE_INI`
export FILE_FASTQ=`ls $SOURCE_FASTQ`
export EXP_ID=`awk '/^exp_id/{print $2}' $FILE_INI`
export FILE_REF=`awk '/^ref_path/{print $2}' $FILE_INI`
export DIR_OUT=$DIR_WORKSPACE/$EXP_ID
export FILE_OUT=$DIR_OUT/$EXP_ID
export FILE_PRIMERFA=${FILE_OUT}_primer.fa
export DIR_LOG=$DIR_OUT/log

# Setup for submission variables
QSUBVARS="-V -e $DIR_LOG -o $DIR_LOG"
HOLD_ID=-1

# Initialization
module load python
mkdir -p $DIR_WORKSPACE $DIR_LOG

# Local: Make the primer fa
#python $MC4CTOOL makeprimerfa $FILE_INI $FILE_PRIMERFA

# Local: Determine amount of reads in the original fastq file
FASTQ_NL=`wc -l $FILE_FASTQ | tail -n 1 |  awk '{print $1}'`
READSPERFILE=20
export LINESPERFILE=$(($READSPERFILE*4))
export NUM_TASKS=$((($FASTQ_NL+$LINESPERFILE-1)/$LINESPERFILE))

HOLD_ID_LIST=""

#echo $(($FASTQWCL/4)) Reads found in $FILE_FASTQ
#echo Becomes $SPLITFILESNUM files

# Single: Split fastq data over several even sized files,
# 	rename reads to unique numbers
HOLD_ID=`qsub $QSUBVARS $DIR_MC4C/sge_scripts/splitfq.sh | awk '{print $3}'`
HOLD_ID_LIST=$HOLD_ID_LIST"splitfq:$HOLD_ID;"

# Array: Map primers to reads
#HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID -t 1:$NUM_TASKS $DIR_MC4C/sge_scripts/bowtie.sh | awk '{print $3}' | cut -f1 -d.`
#HOLD_ID_LIST=$HOLD_ID_LIST"bowtie:$HOLD_ID;"

# Array: Split reads by primers
#HOLD_ID=`qsub $QSUBVARS -hold_jid_ad $HOLD_ID -t 1:$NUM_TASKS $DIR_MC4C/sge_scripts/splitpr.sh | awk '{print $3}' | cut -f1 -d.`
#HOLD_ID_LIST=$HOLD_ID_LIST"splitpr:$HOLD_ID;"

# Array: Split reads by restriction sites
export FILE_SOURCE=block
HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID -t 1:$NUM_TASKS $DIR_MC4C/sge_scripts/splitre.sh | awk '{print $3}' | cut -f1 -d.`
HOLD_ID_LIST=$HOLD_ID_LIST"splitre:$HOLD_ID;"

# Single: Merge split reads to one zipped file
HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID $DIR_MC4C/sge_scripts/mergere.sh | awk '{print $3}'`
HOLD_ID_LIST=$HOLD_ID_LIST"mergere:$HOLD_ID;"

# Single: Map split reads to reference genome
HOLD_ID=`qsub $QSUBVARS -hold_jid $HOLD_ID $DIR_MC4C/sge_scripts/bwa.sh | awk '{print $3}'`
HOLD_ID_LIST=$HOLD_ID_LIST"bwa:$HOLD_ID;"

echo Jobs submitted: $HOLD_ID_LIST
