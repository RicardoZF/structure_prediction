#!/usr/bin/env python
import sys
from pathlib import Path
import string

inp_fasta = sys.argv[1]
out_fasta = sys.argv[2]

nn = 0
# make fasta
with open(inp_fasta) as rf, open(out_fasta, 'w') as wf:
    for line in rf:
        if ">" in line:
            wio =  f">protein|{string.ascii_uppercase[nn]}\n"
            nn += 1
        else:
            wio = line
        wf.write(wio)
