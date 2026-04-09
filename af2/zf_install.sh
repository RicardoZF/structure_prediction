conda create -p ./envs/af2 python=3.11
conda activate ./envs/af2
# conda install -y -c nvidia cuda=12.2.2
conda install -c conda-forge openmm cuda-version=12
conda install -y -c conda-forge pdbfixer
conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2
wget -q -P ./alphafold/common/ \
https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
pip3 install -r requirements.txt --no-cache-dir
pip3 install  jax==0.4.26 jaxlib==0.4.26+cuda12.cudnn89 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
conda install -c conda-forge cudnn=8.9
