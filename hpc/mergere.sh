#$ -l h_rt=00:30:00
#$ -l h_vmem=1G

FILELIST=""
for INDEX in `seq $SPLITFILESNUM`;
do
	FILELIST="$FILELIST ${FILE_OUT}_${INDEX}.splitre.fq"
done

cat ${FILELIST} | gzip > ${FILE_OUT}.splitre.fq.gz
