#$ -cwd
#$ -o logs/run.log
#$ -e logs/run.log
#$ -N run_wg_blimp
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

conda activate ./wg-blimp-env
wg-blimp run-snakemake-from-config \
    --cores 10 \
    --nodes 1 \
    --cluster "qsub -pe local 10 -l mem_free=8G,h_vmem=8G,h_fsize=100G" \
    wg-blimp-config.yaml

echo "**** Job ends ****"
date
