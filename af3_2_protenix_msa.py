#!/usr/bin/env python
import sys
import json
import string
import numpy as np
from pathlib import Path

#if len(sys.argv) < 4:
#    print('usage: ./onvert_a2p_msa.py input_fasta af3_out_json output_dir/')
#    sys.exit()

if len(sys.argv) < 3:
    print('usage: ./onvert_a2p_msa.py af3_out_json output_dir/')
    sys.exit()

#fasta = sys.argv[1]
af3_out_data = sys.argv[1]
output_dir = sys.argv[2]

assert Path(af3_out_data).exists(), f'{af3_out_data} not exists !!!'

output_dir = Path(output_dir).resolve()

msa_dir = output_dir / 'msa'
msa_dir.mkdir(parents=True, exist_ok=True)

with open(af3_out_data) as  f:
    datas = json.load(f)
    seq_obj = datas['sequences']
    
msa_dir_list = []
nn = 1
for so in seq_obj:
    msa_dir_nn = msa_dir / str(nn)
    msa_dir_nn.mkdir(parents=True, exist_ok=True)

    non_pair_file = msa_dir_nn.joinpath('non_pairing.a3m')
    non_pair_data = so['protein']['unpairedMsa']
    non_pair_file.write_text(non_pair_data)

    pair_file = msa_dir_nn.joinpath('pairing.a3m')
    pair_data = so['protein']['pairedMsa']
    pair_file.write_text(pair_data)
    msa_dir_list.append(msa_dir_nn)
    nn += 1
    
#name = Path(fasta).stem
#protenix_inp = output_dir / f'{name}.fasta'
#
#i = 0
#with open(fasta) as rf, open(protenix_inp, 'w') as wf:
#    for line in rf:
#        if ">" in line:
#            msa_dir_nn =  msa_dir_list[i]
#            wio =  f">protein|{string.ascii_uppercase[i]}|{msa_dir_nn}\n"
#            i += 1
#        else:
#            wio = line
#        wf.write(wio)
