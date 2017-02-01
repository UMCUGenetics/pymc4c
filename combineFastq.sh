#IDMATCH=elsewhere
#INFASTQS='./test.FASTQ ./test.FASTQ'
#OUTFASTQ=test.FASTQ.gz

INFASTQS=${@:2}
OUTFASTQ=$1.fastq.gz
IDMATCH=$1.idmatch

awk \
	-v IDMATCH="$IDMATCH" \
	'FNR == 1 { FILENUM += 1; READNUM = 0 }
	!((NR) % 4 == 1) {
		print}
	((NR) % 4 == 1)	{
		READID="@RD:"READNUM";IN:"FILENUM
	        print READID
	        print READID"\tFILE:"FILENAME"\t"$0>IDMATCH;
	        READNUM += 1}' \
	                $INFASTQS | \
	                gzip > $OUTFASTQ
