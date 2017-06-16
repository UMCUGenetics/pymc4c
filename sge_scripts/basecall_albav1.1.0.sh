#$ -l h_rt=05:00:00
#$ -l h_vmem=6G
#$ -pe threaded 6
#$ -cwd
$DEBUG_MODE

cd $SGE_O_WORKDIR
BASECALL=/hpc/cog_bioinf/kloosterman/tools/albacore_v1.1.0/env/bin/read_fast5_basecaller.py
FILE_INPUT=`realpath $1`
FILE_OUT=`realpath $2`
F5_CFG=$3

FOLDER_NAME=$((SGE_TASK_ID-1))

mkdir -p ${FILE_OUT}/${FOLDER_NAME}
echo "Process begins at: `date`" 
echo "Reading from: ${FILE_INPUT}/${FOLDER_NAME}"
echo "Writing to: ${FILE_OUT}/${FOLDER_NAME}"
module load python/3.4.3

$BASECALL \
		--output_format fastq
        --input ${FILE_INPUT}/${FOLDER_NAME}/ \
        --worker_threads 6 \
        --save_path ${FILE_OUT}/${FOLDER_NAME}/ \
        --config ${F5_CFG}

echo "Process finished at: `date`" 