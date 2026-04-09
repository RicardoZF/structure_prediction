#!/sugon_store/pub_data/tools/chai-lab/envs/chai1/bin/python

import os, sys
import numpy as np


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: get_plddt.py output_dpath")
        sys.exit(0)

    inp = sys.argv[1]
    out = os.path.join(inp, "plddt.dat")
    plddt_scores = []

    for i in range(5):
        fpath = os.path.join(inp, f"ranked_{i}.pdb")
        with open(fpath) as lines:
            plddt = [float(x.split()[-2])/100.0 for x in lines \
                    if x.startswith("ATOM") and x.split()[2] == "CA"]
            avg_plddt = np.mean(plddt)
            plddt_scores.append(avg_plddt)

    with open(out, 'w') as tf:
        for i, p in enumerate(plddt_scores):
            tf.write(f'{i} {p:.6f}\n')
