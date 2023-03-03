#$ -cwd
#$ -o logs/run.log
#$ -e logs/run.log
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load nextflow/22.10.7
module load singularity/3.7.0

base_dir=$MYSCRATCH/BiocMAP_benchmarking/methylseq
mkdir -p $base_dir

nextflow run nf-core/methylseq \
    --input $(git rev-parse --show-toplevel)/BiocMAP_benchmark/sample_sheet.csv \
    --outdir $base_dir/out \
    -w $base_dir/work \
    -profile singularity \
    --genome GRCh38

echo "**** Job ends ****"
date
