#IDMATCH=elsewhere
#INFASTAS='./test.fasta ./test.fasta'
#OUTFASTA=test.fasta.gz

INFASTAS=${@:2}
OUTFASTA=$1.fasta.gz
IDMATCH=$1.idmatch

awk \
	-v IDMATCH="$IDMATCH" \
	'FNR == 1 { FILENUM += 1; READNUM = 0 }
	!((NR) % 4 == 1) {
		print}
	((NR) % 4 == 1)	{
		READID=">RD:"READNUM";IN:"FILENUM
	        print READID
	        print READID"\tFILE:"FILENAME"\t"$0>IDMATCH;
	        READNUM += 1}' \
	                $INFASTAS | \
	                gzip > $OUTFASTA
