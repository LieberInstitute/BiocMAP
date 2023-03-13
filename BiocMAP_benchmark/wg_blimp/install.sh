FASTQ_DIR=/dcs04/lieber/lcolladotor/ageNeunSortedWGBS_LIBD001/ageNeunSortedWGBS/raw-data/FASTQ
FASTA=/dcs04/lieber/lcolladotor/BiocMAP_benchmark_LIBD001/BiocMAP_ref_temp/hg38_gencode_v34_main/assembly_hg38_gencode_v34_main.fa
OUT_DIR=/dcs04/lieber/lcolladotor/BiocMAP_benchmark_LIBD001/wg-blimp/out

mkdir -p $OUT_DIR

#   Create a small environment including mamba, since conda is extremely slow
conda create -y -n mamba_env -c conda-forge mamba
conda activate mamba_env

mamba create -p ./wg-blimp-env wg-blimp r-base==4.1.1 mamba
conda activate ./wg-blimp-env

wg-blimp create-config \
    --cores-per-job=4 \
    $FASTQ_DIR/ \
    $FASTA \
    WGC052316L,WGC059613L \
    WGC059614L,WGC052317L \
    $OUT_DIR wg-blimp-config.yaml

#   Then manually modified this attribute in wg-blimp-config.yaml:
# rawsuffixregex:
#   first: _combined_R1\.(fastq|fq)(\.gz)?
#   second: _combined_R2\.(fastq|fq)(\.gz)?
#
#   Also changed aligner to bwa-meth after gemBS had strange issues
