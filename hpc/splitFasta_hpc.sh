#$ -l h_rt=08:00:00
#$ -l h_vmem=8G

#$ -pe threaded 2

#$ -cwd

#$ -o ~/sge
#$ -e ~/sge

#$ -M r.straver-2@umcutrecht.nl
#$ -m bes

BASEDIR=/home/cog/rstraver/Workspace/pymc4c_hpc
EXPNAME='NPC-PCDHa11-NP'

OUTFASTA=$BASEDIR/data/splitbow_${EXPNAME}

#TARGETFA=./data/${EXPNAME}.fa
INFASTAS=/hpc/cog_bioinf/data/aallahyar/My_Works/4C_PacBio/62_Running_ThePipeline_ForAllExperiments/CMB_Files/CMB_NPC-PCDHa11-NP.fasta
#head $INFASTQS > $OUTFASTA.mini.fa
#INFASTQS=$OUTFASTA.mini.fa
#IDMATCH=elsewhere
#INFASTQS='./test.FASTQ ./test.FASTQ'
#OUTFASTQ=test.FASTQ.gz
#65536
#INFASTQS=${@:2}
#OUTFASTA=$1
IDMATCH=$OUTFASTA.idmatch

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
	((NR % 65536) == 1) {
		++C;
		CUROUT = OUTFASTA"_"C".fa" }
	((NR) % 2 == 1)	{
		READID = "RD:"READNUM";IN:"FILENUM;
		print ">"READID > CUROUT;
	        ++READNUM}
	((NR) % 2 == 0) {
		print $0 > CUROUT }
	END { print C }' \
	                $INFASTAS #> \
				#$OUTFASTA.fastq.gz
