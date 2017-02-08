#$ -l h_rt=08:00:00
#$ -l h_vmem=8G

#$ -pe threaded 2

#$ -cwd

#$ -o ~/sge
#$ -e ~/sge

BASEDIR=/home/cog/rstraver/Workspace/pymc4c_hpc
EXPNAME='NPC-PCDHa11-NP'

OUTFASTQ=$BASEDIR/data/splitbowq_${EXPNAME}
INFASTQS=/hpc/cog_bioinf/data/aallahyar/My_Works/4C_PacBio/62_Running_ThePipeline_ForAllExperiments/NPS_Files/NPS_${EXPNAME}_--000009--_Best.fastq
IDMATCH=$OUTFASTQ.idmatch

awk \
	-v IDMATCH="$IDMATCH" \
	-v OUTFASTQ="$OUTFASTQ" \
	'
	BEGIN {
		C=1;
		READNUM = 1;
		CUROUT = OUTFASTQ"_"C".fa";
		ALTOUT = OUTFASTQ"_"C".fq" }
	(FNR == 1) {
		++FILENUM }
	((NR) % 4 == 1)	{
		READID = "RD:"READNUM";IN:"FILENUM;
		print ">"READID > CUROUT
		print ">"READID > ALTOUT;
	        ++READNUM }
	((NR) % 4 == 2) {
		print $0 > CUROUT }
	((NR) % 4 != 1) {
		print $0 > ALTOUT }
	((NR % 65536) == 0) {
		++C;
		CUROUT = OUTFASTQ"_"C".fa"
		ALTOUT = OUTFASTQ"_"C".fq" }
	END { print C }' \
	                $INFASTQS
