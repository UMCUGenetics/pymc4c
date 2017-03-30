#$ -l h_rt=08:00:00
#$ -l h_vmem=8G

#$ -pe threaded 2

FILE_REF=${FILE_OUT}_$SGE_TASK_ID.block.fa
FILE_IDX=${FILE_OUT}_$SGE_TASK_ID.block
FILE_SAM=${FILE_OUT}_$SGE_TASK_ID.block.sam

$BOWTIE-build -f $FILE_REF $FILE_IDX

$BOWTIE --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.50 \
	--rdg 2,1 \
	--rfg 2,1 \
	--mp 3,2 \
	--ma 2 -a -p 2 -f \
	-x $FILE_IDX \
	-U $FILE_PRIMERFA \
	-S $FILE_SAM
