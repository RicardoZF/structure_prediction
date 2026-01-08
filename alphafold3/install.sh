#git clone https://github.com/google-deepmind/alphafold3.git
#cd ./alphafold3

#conda create -n af3 python=3.11
#conda activate af3
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04
conda install -c conda-forge libgcc libgcc-ng  -y
pip3 install -r zf_requirements.txt

pip3 install --no-deps .
