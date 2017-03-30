#$ -l h_rt=00:30:00
#$ -l h_vmem=1G


#BASEDIR=/home/cog/rstraver/Workspace/pymc4c_hpc
#EXPNAME='NPC-PCDHa11-NP'

OUTFASTQ=$FILE_OUT #$BASEDIR/data/splitbowq_${EXPNAME}
#INFASTQS=/hpc/cog_bioinf/data/aallahyar/My_Works/4C_PacBio/62_Running_ThePipeline_ForAllExperiments/NPS_Files/NPS_${EXPNAME}_--000009--_Best.fastq
IDMATCH=$OUTFASTQ.idmatch

awk \
	-v IDMATCH="$IDMATCH" \
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
	END { print C }' \
	                $FILE_FASTQ
