#$ -cwd
#$ -N align_bs
#$ -o logs/align.log
#$ -e logs/align.log
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -t 1-4
#$ -tc 4
#$ -hold_jid build_reference

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

FASTA=$MYSCRATCH/BiocMAP_ref_temp/hg38_gencode_v34_main/assembly_hg38_gencode_v34_main.fa
MANIFEST=$(git rev-parse --show-toplevel)/benchmark/samples.manifest
GENOME_DIR=$(git rev-parse --show-toplevel)/benchmark/bs_seeker3/bs3/reference_genome

#   Get paths to FASTQs for this sample
R1=$(awk "NR == ${SGE_TASK_ID}" $MANIFEST | cut -d $'\t' -f 1)
R2=$(awk "NR == ${SGE_TASK_ID}" $MANIFEST | cut -d $'\t' -f 3)

conda activate $PWD/bs_env
cd bs3
./bs3-align -g $FASTA -1 $R1 -2 $R2 -o test.bam -d $GENOME_DIR

echo "**** Job ends ****"
date
