#$ -cwd
#$ -N build_reference
#$ -o logs/build_reference.log
#$ -e logs/build_reference.log
#$ -l mem_free=40G,h_vmem=40G,h_fsize=100G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

FASTA=$MYSCRATCH/BiocMAP_ref_temp/hg38_gencode_v34_main/assembly_hg38_gencode_v34_main.fa

conda activate $PWD/bs_env
cd bs3
./bs3-build -f $FASTA -L 5

echo "**** Job ends ****"
date
