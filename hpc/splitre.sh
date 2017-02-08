#$ -l h_rt=00:30:00
#$ -l h_vmem=1G

python $MC4CTOOL \
	splitreads \
	$FILE_DATAINFO \
	${FILE_OUT}_$SGE_TASK_ID.splitpr.fq \
	${FILE_OUT}_$SGE_TASK_ID.splitre.fq \
	$EXPNAME
