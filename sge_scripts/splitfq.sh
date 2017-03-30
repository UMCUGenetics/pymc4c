#$ -l h_rt=00:30:00
#$ -l h_vmem=1G

echo "Source: [$FILE_OUT]"
echo "Output: [$FILE_FASTQ]"
echo "Split size: [$LINESPERFILE] lines"


awk \
	-v OUTFASTQ="$FILE_OUT" \
	-v LINESPERFILE="$LINESPERFILE" \
	'
	BEGIN {
		C=1;
		READNUM = 1;
		CUROUT = OUTFASTQ"_"C".block.fa";
		ALTOUT = OUTFASTQ"_"C".block.fq" }
	(FNR == 1) {
		++FILENUM }
	((NR) % 4 == 1)	{
		READID = "Fl.Id:"FILENUM";Rd.Id:"READNUM;
		print ">"READID > CUROUT
		print ">"READID > ALTOUT;
	        ++READNUM }
	((NR) % 4 == 2) {
		print $0 > CUROUT }
	((NR) % 4 != 1) {
		print $0 > ALTOUT }
	((NR % LINESPERFILE) == 0) {
		++C;
		CUROUT = OUTFASTQ"_"C".block.fa"
		ALTOUT = OUTFASTQ"_"C".block.fq" }
	' $FILE_FASTQ
