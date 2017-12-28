#$ -l h_rt=00:30:00
#$ -l h_vmem=1G
$DEBUG_MODE

echo "Source: [$FILE_FASTQ]"
echo "Output: [$FILE_OUT]"
echo "Split size: [$LINESPERFILE] lines"

awk \
	-v OUTFASTQ="$FILE_OUT" \
	-v LINESPERFILE="$LINESPERFILE" \
	'
	BEGIN {
		C=1;
		READNUM = 1;
		CUROUT = OUTFASTQ"_"C".block.fa";
		ALTOUT = OUTFASTQ"_"C".block.fq" 
		READLEN = 0
	}
	(FNR == 1) {
		++FILENUM 
	}
	((FNR) % 4 == 2) {
		++READNUM
		READLEN = length($0)
		READBP = $0
		READID = "Fl.Id:"FILENUM";Rd.Id:"READNUM";Rd.Ln:"length($0)
	}
	((FNR) % 4 == 0 && READLEN >= 500 && READLEN <= 10000) {
		print ">"READID > CUROUT
		print READBP > CUROUT 

		print ">"READID > ALTOUT;
		print READBP > ALTOUT 
		print "+" > ALTOUT
		print $0 > ALTOUT 
	}
	((FNR % LINESPERFILE) == 0) {
		++C;
		CUROUT = OUTFASTQ"_"C".block.fa"
		ALTOUT = OUTFASTQ"_"C".block.fq" 
	}
	' $FILE_FASTQ

# TODO: Add remove read if less than 500 bp