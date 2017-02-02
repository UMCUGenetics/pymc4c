#IDMATCH=elsewhere
#INFASTQS='./test.FASTQ ./test.FASTQ'
#OUTFASTQ=test.FASTQ.gz

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
	((NR % 16384) == 0) {
		++C;
		CUROUT = OUTFASTA"_"C".fa" }
	((NR) % 4 == 1)	{
		READID = "RD:"READNUM";IN:"FILENUM;
		print ">"READID > CUROUT;
	        print ">"READID"\t"FILENAME"\t"$0 | "gzip >" IDMATCH".gz";
		print "@"READID | "gzip >" OUTFASTA".fastq.gz"
	        ++READNUM}
	((NR) % 4 == 2) {
		print $0 > CUROUT }
	!((NR) % 4 == 1) {
		print $0 | "gzip >" OUTFASTA".fastq.gz"}
	END { print C }' \
	                $INFASTQS #> \
				#$OUTFASTA.fastq.gz
