#$ -l h_rt=00:30:00
#$ -l h_vmem=1G
$DEBUG_MODE

python $MC4CTOOL \
	cleavereads \
	$FILE_INI \
	${FILE_OUT}_$SGE_TASK_ID.sam \
	${FILE_OUT}_$SGE_TASK_ID.block.fq \
	${FILE_OUT}_$SGE_TASK_ID.splitpr.fq
