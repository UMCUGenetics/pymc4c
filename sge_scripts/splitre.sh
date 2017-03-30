#$ -l h_rt=00:30:00
#$ -l h_vmem=1G


if [ ! -n "FILE_SOURCE" ]; then
    FILE_SOURCE=splitpr
fi

python $MC4CTOOL \
	splitreads \
	$FILE_INI \
	${FILE_OUT}_$SGE_TASK_ID.$FILE_SOURCE.fq \
	${FILE_OUT}_$SGE_TASK_ID.splitre.fq
