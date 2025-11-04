import sys
import linecache
import numpy as np

fastaf = sys.argv[1]
mut_list = sys.argv[2]

muts = np.loadtxt(mut_list, dtype=str)
seq = ''.join([i.split()[0] for i in linecache.getlines(fastaf)[1:]])

for mut in muts:
  p = eval(mut[1:-1])
  m = mut[-1]
  
  assert seq[p-1] == mut[0]
  
  seq2 = seq[:p-1] + mut + seq[p:]
  
  with open(mut + '.fasta', 'w') as f:
    f.write('>{}\n{}\n'.format(mut, seq2))
