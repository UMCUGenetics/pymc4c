BOWTIE=/hpc/cog_bioinf/data/aallahyar/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2

EXPNAME='FL-HS3'

for PART in `seq 0 3`
do
	TARGETFA=./data/${EXPNAME}_${PART}.fa
	TARGETBT=./data/${EXPNAME}_${PART}
	PRIMERFA=./local/${EXPNAME}.fa
	TARGETSAM=./data/${EXPNAME}_${PART}.sam

	$BOWTIE-build -f $TARGETFA $TARGETBT

	$BOWTIE --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.50 \
		--rdg 2,1 \
		--rfg 2,1 \
		--mp 3,2 \
		--ma 2 -a -p 2 -f \
		-x $TARGETBT \
		-U $PRIMERFA \
		-S $TARGETSAM
done
