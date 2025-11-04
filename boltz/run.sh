#!/bin/bash

inp_fasta=$1
output=$2
module load cuda12.4

# The directory where to download the data and model.
CACHE="/sugon_store/pub_data/databases/Models/boltz"

PKG="$(dirname `realpath $0`)"
boltz=$PKG/src/boltz/main.py
PYTHON=$PKG/envs/boltz_v2/bin/python3.11
export PYTHONPATH=$PKG/src:$PYTHONPATH

inp_file=$output/out.yaml
if [[ "${inp_fasta##*.}" == "yaml" ]];then
  cp $inp_fasta $inp_file
else
  $PYTHON $PKG/make_msa_fasta.py $inp_fasta $output
fi

cmd="$PYTHON $boltz predict $inp_file --out_dir $output --output_format pdb --cache $CACHE --override --diffusion_samples 5"
echo "run boltz: $cmd"
$cmd
