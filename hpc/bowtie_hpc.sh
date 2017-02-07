#$ -l h_rt=08:00:00
#$ -l h_vmem=8G

#$ -pe threaded 2

#$ -cwd

#$ -o ~/sge
#$ -e ~/sge

#$ -M r.straver-2@umcutrecht.nl
#$ -m bes

BOWTIE=/hpc/cog_bioinf/data/aallahyar/My_Works/Useful_Sample_Codes/Bowtie2/bowtie2-2.2.6/bowtie2

BASEDIR=/home/cog/rstraver/Workspace/pymc4c_hpc
EXPNAME='NPC-PCDHa11-NP'


#TARGETFA=./data/${EXPNAME}.fa
TARGETFA=/hpc/cog_bioinf/data/aallahyar/My_Works/4C_PacBio/62_Running_ThePipeline_ForAllExperiments/CMB_Files/CMB_NPC-PCDHa11-NP.fasta
TARGETBT=$BASEDIR/data/bigbow_${EXPNAME}
PRIMERFA=$BASEDIR/local/${EXPNAME}.fa
TARGETSAM=$BASEDIR/data/bigbow_${EXPNAME}.sam

$BOWTIE-build -f $TARGETFA $TARGETBT

$BOWTIE --local -D 20 -R 3 -N 0 -L 15 -i S,1,0.50 \
	--rdg 2,1 \
	--rfg 2,1 \
	--mp 3,2 \
	--ma 2 -a -p 2 -f \
	-x $TARGETBT \
	-U $PRIMERFA \
	-S $TARGETSAM
