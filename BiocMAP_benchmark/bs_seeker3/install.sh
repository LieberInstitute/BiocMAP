#   Create a small environment including mamba, since conda extremely slow
conda create -y -n mamba_env -c conda-forge mamba
conda activate mamba_env

#   Install python dependencies in new conda environment
mamba create -y -p $PWD/bs_env python=2.7
conda activate $PWD/bs_env
pip install pysam Matplotlib

#   Clone main repo
git clone https://github.com/khuang28jhu/bs3

#   Clone SNAP 1.0.0, a supposed dependency
wget https://github.com/amplab/snap/archive/refs/tags/v1.0.0.tar.gz
tar -xzf v1.0.0.tar.gz
rm v1.0.0.tar.gz

#   Test basic functionality
cd bs3
./bs3-build -h

conda deactivate
