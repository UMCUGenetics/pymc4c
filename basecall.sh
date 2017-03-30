#$ -l h_rt=04:00:00
#$ -l h_vmem=24G
#$ -pe threaded 12

BASECALL=/hpc/cog_bioinf/kloosterman/users/wkloosterman/tools/albacore/albacore_venv/bin/read_fast5_basecaller.py
FILE_INPUT=$1
FILE_OUT=${FILE_INPUT}_basecalled
F5_CFG=FLO-MIN106_LSK108_linear.cfg

THIS_ID=$((SGE_TASK_ID-1))

mkdir -p ${FILE_OUT}/${THIS_ID}
echo ${FILE_INPUT}/${THIS_ID}
module load python/3.4.3

$BASECALL \
        --input ${FILE_INPUT}/${THIS_ID}/ \
        --worker_threads 12 \
        --save_path ${FILE_OUT}/${THIS_ID}/ \
        --config ${F5_CFG}