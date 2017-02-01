set -e

export PYTHON=python
export MAIN=mc4c.py
export BOWTIE=bowtie2
export DIR_TOOL=./
export DIR_DATA=./data
export DIR_IN=$DIR_DATA/in
export DIR_RFS=$DIR_DATA/rfs
export DIR_CMB=$DIR_DATA/cmb
export DIR_LOG=$DIR_DATA/log
export FILE_DATAINFO=./Dataset_info.mat

export MAX_CPU=4


#LIST_INDEXES=`$PYTHON mc4c.py listindex $FILE_DATAINFO`

mkdir -p $DIR_RFS $DIR_CMB $DIR_LOG


for INDEX in Upp2 #`cat indexList`
do
	echo Working on $INDEX
	FILE_RFS=$DIR_RFS/$INDEX.rfs.fa
	FILE_CMB=$DIR_CMB/$INDEX.fasta
	FILE_REF=$DIR_RFS/$INDEX.rfs.ref
	FILE_SAM=$DIR_RFS/$INDEX.sam

	echo RFS Fasta creation...
	if [ ! -f $FILE_RFS ];
	then
		$PYTHON $MAIN makeprimerfa \
			$FILE_DATAINFO \
			$FILE_RFS \
			$INDEX \
				> $DIR_LOG/$INDEX.rfsfasta.log
	fi

	echo Bowtie2-build in action...
	if [ ! -f $FILE_REF ];
	then
		echo $BOWTIE-build \
			-f \
			$FILE_CMB \
			$FILE_REF \
				> $DIR_LOG/$INDEX.btbuild.log
		touch $FILE_REF
	fi


	echo Bowtie2 in action...
	if [ ! -f $FILE_SAM ];
	then
		echo $BOWTIE \
			--local \
			-D 20 \
			-R 3 \
			-N 0 \
			-L 15 \
			-i S,1,0.50 \
			--rdg 2,1 \
			--rfg 2,1 \
			--mp 3,2 \
			--ma 2 \
			-a \
			-p $MAX_CPU \
			-f \
			-x $FILE_REF \
			-U $FILE_RFS \
			-S $FILE_SAM \
				> $DIR_LOG/$INDEX.btsam.log

		touch $FILE_SAM
	fi

	#samtools view -S RFS_LVR-HS2-NP2.sam -b > RFS_LVR-HS2-NP2.bam
	#samtools view RFS_LVR-HS2-NP2.sortname.bam
done
