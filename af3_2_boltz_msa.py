#!/usr/bin/env python
import sys
import json
import string
import numpy as np
from pathlib import Path

if len(sys.argv) < 3:
    print('usage: ./af3_2_boltz_msa.py af3_out_json output_dir/')
    sys.exit()

af3_out_data = sys.argv[1]
output_dir = sys.argv[2]

assert Path(af3_out_data).exists(), f'{af3_out_data} not exists !!!'

output_dir = Path(output_dir).resolve()
output_dir.mkdir(exist_ok=True)

with open(af3_out_data) as  f:
    datas = json.load(f)
    seq_obj = datas['sequences']
    
msa_dir_list = []
nn = 0
for so in seq_obj:
    if 'protein' in so:
        non_pair_data = so['protein']['unpairedMsa']
        pair_data = so['protein']['pairedMsa']
        msa_file = output_dir.joinpath(f'{nn}.a3m')
        msa_file.write_text(non_pair_data+pair_data)
        nn += 1
