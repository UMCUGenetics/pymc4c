#$ -l h_rt=00:30:00
#$ -l h_vmem=1G
#$ -pe threaded 12

BASECALL=/hpc/cog_bioinf/kloosterman/users/wkloosterman/tools/albacore/albacore_venv/bin/read_fast5_basecaller.py
FILE_INPUT=/hpc/cog_bioinf/kloosterman/users/mroosmalen/scripts/NAP/testing/input
FILE_OUT=/hpc/cog_bioinf/data/rstraver/basecall/
F5_CFG=FLO-MIN106_LSK108_linear.cfg

mkdir -p ${FILE_OUT}/$SGE_TASK_ID

python $BASECALL \
	--input ${FILE_INPUT}/$SGE_TASK_ID \
	--worker_threads 12 \
	--save_path ${FILE_OUT}/$SGE_TASK_ID/ \
	--config ${F5_CFG}
