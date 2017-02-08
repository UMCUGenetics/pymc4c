set -e

export DIR_TOOL="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PYTHON=python
export MAIN=mc4c.py
export BOWTIE=bowtie2
export DIR_DATA=./data
export DIR_IN=$DIR_DATA/in
export DIR_RFS=$DIR_DATA/rfs
export DIR_CMB=$DIR_DATA/cmb
export DIR_LOG=$DIR_DATA/log
export FILE_DATAINFO=./Dataset_info.mat

export MAX_CPU=4

QSUBVARS="-V -e $DIR_LOG -o $DIR_LOG"

#LIST_INDEXES=`$PYTHON mc4c.py listindex $FILE_DATAINFO`

mkdir -p $DIR_RFS $DIR_CMB $DIR_LOG

module load python

# Split fastq data over several even sized files - Single
# Rename read names to unique numbers
echo qsub $QSUBVARS splitFastq_hpc.sh

# Map primers to reads - Array
echo qsub $QSUBVARS bowtie_hpc_array.sh

# Split reads by primers - Array
echo qsub $QSUBVARS cleavereads.sh
#python mc4c.py cleavereads ./local/Dataset_info.mat ./local/splitbowq_NPC-PCDHa11-NP_32.sam ./local/splitbowq_NPC-PCDHa11-NP_32.fq ./local/splitbowq_NPC-PCDHa11-NP_32.prcleave.fq NPC-PCDHa11-NP

# Split reads by restriction sites - Array
echo qsub $QSUBVARS splitreads.sh
#python mc4c.py splitreads ./local/Dataset_info.mat ./local/splitbowq_NPC-PCDHa11-NP_32.prcleave.fq ./local/splitbowq_NPC-PCDHa11-NP_32.split.fq NPC-PCDHa11-NP

# Merge split reads - Single
echo qsub $QSUBVARS "Merge step"

# Map split reads to reference genome - Single
echo qsub $QSUBVARS "Map step"
