#$ -l h_rt=00:30:00
#$ -l h_vmem=1G

cat ${FILE_OUT}_*.splitre.fq | gzip > ${FILE_OUT}.splitre.fq.gz
