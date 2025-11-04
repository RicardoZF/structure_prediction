#!/bin/bash

input_file=$1
output=$2

if [[ -f $output/out/out_model.cif ]];then
  echo "output:$output/out/out_model.cif is exists, exit !!!"
  exit
fi

db_dir=/sugon_store/pub_data/databases/alphafold3
model_dir=/sugon_store/pub_data/databases/Models/alphafold3

PKG="/sugon_store/pub_data/tools/alphafold3"
#PKG="$(dirname `realpath $0`)"
#echo "PKG: $PKG"
AF3_ENV="$PKG/envs/af3"
PYTHON_EXE=$AF3_ENV/bin/python3.11
export PATH=$AF3_ENV/bin:$PATH
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter" 


# ------------------  解决显存溢出问题 ------------------ #
# 统一内存: 启用统一内存允许程序在没有足够空间时将GPU内存溢出到主机内存。这可以防止OOM，但代价是程序变慢，因为访问主机内存而不是设备内存。

# 这会禁用预分配行为。JAX 将改为根据需要分配 GPU 内存，从而可能减少总体内存使用量。但是，此行为更容易导致 GPU 内存碎片，这意味着使用大部分可用 GPU 内存的 JAX 程序在禁用预分配的情况下可能会出现 OOM。
export XLA_PYTHON_CLIENT_PREALLOCATE=false
# 允许 CPU 和 GPU 共享内存空间，避免 GPU OOM（Out of Memory）错误
export TF_FORCE_UNIFIED_MEMORY=true
# 如果启用了预分配，这将使 JAX 预分配 XX% 的 GPU 总内存，而不是默认的 75%。降低预分配量可以修复 JAX 程序启动时发生的 OOM。
export XLA_CLIENT_MEM_FRACTION=3.2


mkdir -pv $output
output="`realpath $output`"
input_json_path="$output/inp.json"
if [[ "${input_file##*.}" == "json" ]];then
  cp $input_file $input_json_path
else
   # Convert fasta to Protenix json
  $PYTHON_EXE $PKG/fasta.py $input_file $output
fi

$PYTHON_EXE $PKG/run_alphafold.py --model_dir=$model_dir --db_dir=$db_dir --flash_attention_implementation=xla --json_path=$input_json_path --output_dir=$output
