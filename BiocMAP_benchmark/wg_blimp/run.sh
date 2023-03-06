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
    --cores 2 \
    --nodes 10 \
    --cluster "qsub -pe local 2 -l mem_free=20G,h_vmem=20G,h_fsize=100G" \
    wg-blimp-config.yaml

echo "**** Job ends ****"
date

#   Tried but got strange errors:
# --cores 2 \
# --nodes 10 \
# --cluster "qsub -pe local 2 -l mem_free=20G,h_vmem=20G,h_fsize=100G" \
