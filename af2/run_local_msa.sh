#!/bin/bash

if [ $# -ne 2 ]; then   
  echo "Usage: run_local_msa.sh input.fasta output/"
  exit 0;
fi 

inp=`realpath $1` 
out=$2

if [ ! -f $inp ]; then
  echo "[INFO]`date`: input fasta file $inp not exists, exit now!!!" 
  exit 0; 
fi 

# Python and binary environment
#PKG_ROOT=`realpath $0`
export PKG_ROOT=/sugon_store/pub_data/tools/alphafold
echo "[INFO]`date`: package root directory $PKG_ROOT"
export AF2_ENV="$PKG_ROOT/envs/af2"
export AF2_DATA="/sugon_store/pub_data/databases/alphafold3"
export PATH="$PKG_ROOT/scripts:$AF2_ENV/bin:$PATH"
export PYTHON_EXE="$AF2_ENV/bin/python" #`which python`
export SEED="1"
echo "[INFO]`date`: python environment $PYTHON_EXE"

# ----------------------- Define a MSA calculate function -----------------------
# MSA for single chain
function generate_msa(){
  # chainA.fasta
  inp_fas=$1 
  # output/xxx/msas/A
  out_dir=$2 
  # make a temp file token 
  tmptoken=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c10)
  tmpoutdir=$3/msa_$tmptoken
  #tmpoutdir="$out_dir/tmp/msa_$tmptoken"
  mkdir -pv $tmpoutdir

  # UNIREF90
  if [ ! -f $out_dir/uniref90_hits.sto ]; then 
    cmd="jackhmmer -o /dev/null -A $tmpoutdir/uniref90_hits.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 $inp_fas $AF2_DATA/uniref90.fasta" 
    echo "[INFO]`date`: running cmd: $cmd"
    $cmd & 
  fi

  # MGY
  if [ ! -f $out_dir/mgnify_hits.sto ]; then
    cmd="jackhmmer -o /dev/null -A $tmpoutdir/mgnify_hits.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 $inp_fas $AF2_DATA/mgy_clusters.fa"
    echo "[INFO]`date`: running cmd: $cmd"
    $cmd &
  fi

  # BFD
  if [ ! -f $out_dir/bfd_uniclust_hits.a3m ]; then
    cmd="hhblits -i $inp_fas -cpu 8 -oa3m $tmpoutdir/bfd_uniclust_hits.a3m -n 3 -e 0.001 -maxseq 1000000 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -d $AF2_DATA/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt -d $AF2_DATA/uniref30/UniRef30_2021_03"
    echo "[INFO]`date`: running cmd: $cmd"
    $cmd & 
  fi

  sleep 120s
  # UNIPROT
  if [ ! -f $out_dir/uniprot_hits.sto ]; then
    cmd="jackhmmer -o /dev/null -A $tmpoutdir/uniprot_hits.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 $inp_fas $AF2_DATA/uniprot.fasta"
    echo "[INFO]`date`: running cmd: $cmd"
    $cmd & 
  fi
  wait 
  
  # now processing templates 
  if [ ! -f $out_dir/pdb_hits.sto ]; then
     if [ ! -f $tmpoutdir/uniref90_hits.sto ]; then 
       cp -uv $out_dir/uniref90_hits.sto $tmpoutdir/uniref90_hits.sto
     fi 

     cmd="hmmbuild  --cpu 8 --hand --amino $tmpoutdir/uniref90.hmmer $tmpoutdir/uniref90_hits.sto"
     $cmd & 
     echo "[INFO]`date`: running cmd: $cmd"
     wait
     cmd="hmmsearch --noali --seed $SEED --cpu 8 --F1 0.1 --F2 0.1 --F3 0.1 --incE 100 -E 100 --domE 100 --incdomE 100 -A $tmpoutdir/pdb_hits.sto $tmpoutdir/uniref90.hmmer $AF2_DATA/pdb_seqres/pdb_seqres.txt"
     echo "[INFO]`date`: running cmd: $cmd"
     $cmd & 
     wait
     rm -rf $tmpoutdir/uniref90.hmmer
  fi  
  
  cp -uv $tmpoutdir/* $out_dir
  echo "[INFO]`date`: Complete MSA searching for $inp_fas $out_dir"
  rm -rfv $tmpoutdir
}

# ----------------------- Define a MSA calculate function -----------------------
for item in `$PYTHON_EXE $PKG_ROOT/scripts/process_fasta.py $inp`
do
  # start MSA
  chainid=`echo $item | cut -d "|" -f 1`
  seq=`echo $item | cut -d "|" -f 2`
  tmptoken=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c10)
  mkdir -pv $out/tmp
  fasta="$out/tmp/chain${chainid}.fasta"
  if [[ ! -f $fasta || ! -s $fasta ]];then
    echo ">chain${chainid}"  >> $fasta
    echo "$seq" >> $fasta
  fi

  mkdir -pv $out/msas/$chainid
  for f in `ls $out/msas/$chainid`
  do
    lc=`cat $out/msas/$chainid/$f | wc -l`
    if [ $lc -lt 1 ]; then
      echo "[INFO] `date`: empty file $out/msas/$chainid/$f, remove it!!!"
      rm -v $out/msas/$chainid/$f
    fi
  done
  generate_msa $fasta $out/msas/$chainid $out/tmp
  echo "[INFO]`date`: finish MSA for chain $chainid !!!"
done
echo "[INFO] `date`: Complete Whole MSA generation process!!!"

exit 999;
