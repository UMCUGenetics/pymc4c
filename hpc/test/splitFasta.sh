#IDMATCH=elsewhere
#INFASTQS='./test.FASTQ ./test.FASTQ'
#OUTFASTQ=test.FASTQ.gz
#65536
INFASTQS=${@:2}
OUTFASTA=$1
IDMATCH=$1.idmatch

awk \
	-v IDMATCH="$IDMATCH" \
	-v OUTFASTA="$OUTFASTA" \
	'
	BEGIN{
		C=0;
		CUROUT = OUTFASTA"_"C".fa"}
	(FNR == 1) {
		++FILENUM;
		READNUM = 0 }
	((NR % 4) == 1) {
		++C;
		CUROUT = OUTFASTA"_"C".fa" }
	((NR) % 2 == 1)	{
		READID = "RD:"READNUM";IN:"FILENUM;
		print ">"READID > CUROUT;
	        ++READNUM}
	((NR) % 2 == 0) {
		print $0 > CUROUT }
	END { print C }' \
	                $INFASTQS #> \
				#$OUTFASTA.fastq.gz
