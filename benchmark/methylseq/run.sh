#$ -cwd
#$ -o logs/run_2.log
#$ -e logs/run_2.log
#$ -l mem_free=15G,h_vmem=15G,h_fsize=1000G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load nextflow/22.10.7
module load singularity/3.6.0

base_dir=/dcs04/lieber/lcolladotor/BiocMAP_benchmark_LIBD001/methylseq
mkdir -p $base_dir

nextflow run nf-core/methylseq \
    --input $(git rev-parse --show-toplevel)/BiocMAP_benchmark/sample_sheet.csv \
    --outdir $base_dir/out \
    -w $base_dir/work \
    -profile singularity \
    --genome GRCh38 \
    -with-report "logs/exec_report_2.html" \
    -resume

echo "**** Job ends ****"
date
