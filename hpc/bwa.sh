#$ -l h_rt=08:00:00
#$ -l h_vmem=32G

#$ -pe threaded 12

$BWA bwasw \
	-b 5 \
	-q 2 \
	-r 1 \
	-z 10 \
	-T 15 \
	-t 12 \
 	-f ${FILE_OUT}.bwa.sam \
	$FILE_REF \
	${FILE_OUT}.splitre.fq.gz
