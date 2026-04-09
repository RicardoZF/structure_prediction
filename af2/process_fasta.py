#!/usr/bin/env python

import os, sys

def read_fasta(fpath):
    seq_list = []
    with open(fpath) as lines:
        seq = ""
        for l in lines:
            if ">" in l:
                if seq != "":
                    seq_list.append(seq)
                seq = ""
            else:
                seq += l.strip("\n").strip() 
        seq_list.append(seq)

    return seq_list


if __name__ == "__main__":
    inp = sys.argv[1]

    chains = "ABCDEFGHIGKLMNOPQRSTUVWXYZ"
    seq_list = read_fasta(inp)

    mapping_dict = {}
    c = 0
    for i, seq in enumerate(seq_list):
        if seq not in mapping_dict:
            mapping_dict[seq] = chains[c]
            c += 1
            
        print("{}|{}".format(mapping_dict[seq], seq))
            
