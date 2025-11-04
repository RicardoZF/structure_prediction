# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Attribution-NonCommercial 4.0 International
# License (the "License"); you may not use this file except in
# compliance with the License. You may obtain a copy of the
# License at

#     https://creativecommons.org/licenses/by-nc/4.0/

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# zf load cuda
#module load gcc/11.2.0
module load cuda12.4
#export CUDA_HOME=/opt/ohpc/pub/apps/cuda12.4

PKG=$(dirname `realpath $0`)
protenix_env="$PKG/envs/protenix"
#protenix_env="/sugon_store/zhoufan/.miniconda3/envs/protenix2"

# Protenix DataBases
export PROTENIX_DATA_ROOT_DIR="/sugon_store/pub_data/databases/ProtenixData"

export LAYERNORM_TYPE=fast_layernorm
# zf add 
export CUTLASS_PATH=$PKG/cutlass
export PATH=$protenix_env/bin:$PATH
PYTHON_EXE=$protenix_env/bin/python3


# Convert fasta to Protenix json
input_file=$1
output=$2
if [[ -f $output/out/seed_101/predictions/out_seed_101_sample_0.cif ]];then
   echo "output: $output/out/seed_101/predictions/out_seed_101_sample_0.cif is exists, exit !!!"
   exit
fi

mkdir -pv $output
input_json_path="$output/inp.json"
if [[ "${input_file##*.}" == "json" ]];then
  cp $input_file $input_json_path
else
  $PYTHON_EXE $PKG/fasta.py $input_file $output
fi

N_sample=5
N_step=200
N_cycle=10
seed=101
#use_deepspeed_evo_attention=true
use_deepspeed_evo_attention=false
# wget -P https://af3-dev.tos-cn-beijing.volces.com/release_model/model_v0.2.0.pt 
#load_checkpoint_path="/af3-dev/release_model/model_v0.2.0.pt"
load_checkpoint_path="/sugon_store/pub_data/databases/ProtenixData/model_v0.5.0.pt"

$PYTHON_EXE $PKG/runner/inference.py \
--seeds ${seed} \
--load_checkpoint_path ${load_checkpoint_path} \
--dump_dir ${output} \
--input_json_path ${input_json_path} \
--use_deepspeed_evo_attention ${use_deepspeed_evo_attention} \
--model.N_cycle ${N_cycle} \
--sample_diffusion.N_sample ${N_sample} \
--sample_diffusion.N_step ${N_step}
