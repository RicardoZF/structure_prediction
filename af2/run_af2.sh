#!/bin/bash

if [ $# -ne 2 ]; then   
  echo "Usage: run_af2multimer.sh input.fasta output/"
  exit 0;
fi  

fasta=`realpath $1` 
if [ ! -f $fasta ]; then
  echo "[INFO]`date`: input fasta file $fasta not exists, exit now!!!" 
  exit 0; 
fi 

# define output directory 
mkdir -pv $2
outdir=`realpath $2`

# CUDA environment
#module load cudnn/7.6.5+cuda10.1 
# no need to set

# Python and binary environment
PKG_ROOT=`realpath $0`
export PKG_ROOT=`dirname $PKG_ROOT`
echo "[INFO]`date`: package root directory $PKG_ROOT"
export AF2_ENV="$PKG_ROOT/envs/af2"
export AF2_DATA="/sugon_store/pub_data/databases/alphafold3"
export PATH="$PKG_ROOT/scripts:$AF2_ENV/bin:$PATH"
export PYTHON_EXE="$AF2_ENV/bin/python" #`which python`
export PYTHONPATH=$AF2_ENV/lib/python3.11/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=$AF2_ENV/lib:$LD_LIBRARY_PATH
echo "[INFO]`date`: python environment $PYTHON_EXE"
# E2E and AmberRelaxation
#---- from alphafold2 run_docker.py -----#
# The following flags allow us to make predictions on proteins that
# would typically be too long to fit into GPU memory.
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
echo "[INFO]`date`: starting AF2 E2E step ..."

# try 3 times 
ntrial=0

while [ $ntrial -lt 3 ]; do
  if [ -f "$2/`basename $fasta | cut -d "." -f 1`/ranked_0.pdb" ]; then
    echo "find previous output $2/`basename $fasta | cut -d "." -f 1`/ranked_0.pdb, exit now!"
    exit 999;
  fi
  
  seqn=`grep ">" $fasta | wc -l`
  if [[ $seqn -eq 1 ]];then
    model_preset="monomer"
    # model_str="model_1,model_2,model_3,model_4,model_5"
  else
    model_preset="multimer"
    # model_str="model_1_multimer,model_2_multimer,model_3_multimer,model_4_multimer,model_5_multimer"
  fi

  $PYTHON_EXE $PKG_ROOT/run_alphafold.py \
     --fasta_paths="$fasta" \
     --output_dir="$outdir" \
     --data_dir="$AF2_DATA" \
     --uniref90_database_path="$AF2_DATA/uniref90.fasta" \
     --mgnify_database_path="$AF2_DATA/mgy_clusters.fa" \
     --template_mmcif_dir="$AF2_DATA/mmcif_files" \
     --max_template_date="2019-12-31" \
     --obsolete_pdbs_path="$AF2_DATA/obsolete.dat" \
     --bfd_database_path="$AF2_DATA/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt" \
     --uniref30_database_path="$AF2_DATA/uniref30/UniRef30_2021_03" \
     --model_preset="$model_preset" \
     --pdb_seqres_database_path="$AF2_DATA/pdb_seqres/pdb_seqres.txt" \
     --uniprot_database_path="$AF2_DATA/uniprot.fasta" \
     --pdb70_database_path="$AF2_DATA/pdb70/pdb70" \
     --use_gpu_relax="True" \
     --num_multimer_predictions_per_model=1
     #--model_names=$model_str \
  ntrial=$(expr $ntrial + 1)
  echo "[INFO] `date`: process the E2E ntrial $ntrial"
done
exit 999;
