#!/usr/bin/env python
import sys
import os
from pathlib import Path
import string
import numpy as np
import yaml

from collections.abc import Mapping
from pathlib import Path

from Bio import SeqIO
from rdkit.Chem.rdchem import Mol

ALPHABET_LIST = string.ascii_uppercase

# define modifications ccd code
replacements = {
    '(HYP)': 'P',
    '(SME)': 'M',
    '(LYZ)': 'K'
}

def replace_and_track(input_seq):
    result = []
    positions = []
    i = 0
    while i < len(input_seq):
        replaced = False
        for tag, replacement in replacements.items():
            if input_seq.startswith(tag, i):
                result.append(replacement)
                positions.append((tag.strip('()'), len(''.join(result))))
                i += len(tag)
                replaced = True
                break
        if not replaced:
            result.append(input_seq[i])
            i += 1
    return ''.join(result), positions


# query msa
def query_msa(fasta, output_dir):
    with Path(fasta).open("r") as f:
        records = list(SeqIO.parse(f, "fasta"))
        
    query_fasta_io = []
    a3m_map = {}
    i = 0
    for seq_record in records:
        header = seq_record.id.split("|")
        seq = str(seq_record.seq)
        if header[0].lower() == "protein" and len(header) < 3:
            seq = replace_and_track(seq)[0]
            out_a3m = Path(output_dir) / f'{i}.a3m'
            a3m_map[seq] = str(out_a3m)
            if not out_a3m.exists():
                query_fasta_io.append(f">query\n{seq}\n")
            i += 1

    if query_fasta_io:
        query_fasta = os.path.join(output_dir, 'query.fasta')
        with open(query_fasta, 'w') as f:
            f.writelines(query_fasta_io)

        a3mpy = str(Path(__file__).parent.resolve() / 'af3_jackhmmer.py')
        os.system(f'{a3mpy} --fasta_path {query_fasta} --output_dir {output_dir}')
    return a3m_map


# fasta to yaml
def parse_fasta(fasta_file, output_dir, a3m_map = None):
    """Parse a fasta file.

    The name of the fasta file is used as the name of this job.
    We rely on the fasta record id to determine the entity type.

    > CHAIN_ID|ENTITY_TYPE
    SEQUENCE
    > CHAIN_ID|ENTITY_TYPE
    ...

    Where ENTITY_TYPE is either protein, rna, dna, ccd or smiles,
    and CHAIN_ID is the chain identifier, which should be unique.
    The MSA_DIR is optional and should only be used on proteins.

    Parameters
    ----------
    fasta_file : Path to the fasta file.
    output_dir: Directory to save output yaml file.

    """
    
    # Read fasta file
    with Path(fasta_file).open("r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    # Make sure all records have a chain id and entity
    for seq_record in records:
        if "|" not in seq_record.id:
            msg = f"Invalid record id: {seq_record.id}"
            raise ValueError(msg)

        header = seq_record.id.split("|")
        assert len(header) >= 2, f"Invalid record id: {seq_record.id}"

        # zf update  
        #chain_id, entity_type = header[:2]
        entity_type, chain_id= header[:2]
        if entity_type.lower() not in {"protein", "dna", "rna", "ligand"}:
            msg = f"Invalid entity type: {entity_type}"
            raise ValueError(msg)
        if chain_id == "":
            msg = "Empty chain id in input fasta!"
            raise ValueError(msg)
        if entity_type == "":
            msg = "Empty entity type in input fasta!"
            raise ValueError(msg)

    # Convert to yaml format
    sequences = []
    chain_ndx = 0
    binder = None
    for seq_record in records:
        # Get chain id, entity type and sequence
        header = seq_record.id.split("|")

        # zf update  
        #chain_id, entity_type = header[:2]
        entity_type, chain_id= header[:2]
        if chain_id not in ALPHABET_LIST: 
            chain_id = ALPHABET_LIST[chain_ndx]
        chain_ndx += 1

        entity_type = entity_type.upper()
        seq = str(seq_record.seq)

        if entity_type == "PROTEIN":            
            seq, ptm_datas = replace_and_track(seq)
            modifications = [{"ccd": ptm[0], "position": ptm[1]} for  ptm in ptm_datas]
            
            if len(header) == 3 and header[2] != "":
                msa_id = header[2]
            else:
                msa_id = a3m_map[seq] if a3m_map else None
            
            molecule = {
                "protein": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": modifications,
                    "msa": msa_id,
                },
            }
        elif entity_type == "RNA":
            molecule = {
                "rna": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": [],
                },
            }
        elif entity_type == "DNA":
            molecule = {
                "dna": {
                    "id": chain_id,
                    "sequence": seq,
                    "modifications": [],
                }
            }
        elif entity_type.upper() == "LIGAND":                
            molecule = {
                "ligand": {
                    "id": chain_id,
                }
            }
            if len(seq) == 3:
                molecule['ligand']["ccd"] = seq
            else:
                molecule['ligand']["smiles"] = seq
                        
            if len(header) == 3 and header[2] == "affinity":
                binder = chain_id

        #elif entity_type.upper() == "CCD":                
        #    molecule = {
        #        "ligand": {
        #            "id": chain_id,
        #            "ccd": seq,
        #        }
        #    }
        #                
        #    if len(header) == 3 and header[2] == "affinity":
        #        binder = chain_id
        #    
        #elif entity_type.upper() == "SMILES":
        #    molecule = {
        #        "ligand": {
        #            "id": chain_id,
        #            "smiles": seq,
        #        }
        #    }
        #                
        #    if len(header) == 3 and header[2] == "affinity":
        #        binder = chain_id

        sequences.append(molecule)

    data = {
        "sequences": sequences,
        "version": 1,
    }
    
     # Only one single molecule can be specified for affinity computation
    if binder:
        data["properties"] = [{"affinity":{"binder": binder}}]

    with open(Path(f'{output_dir}/out.yaml'), 'w') as f:
        yaml.dump(data, f)

if __name__ == '__main__':
    fasta = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.realpath(output_dir)

    a3m_map = query_msa(fasta, output_dir)
        
    parse_fasta(fasta, output_dir, a3m_map = a3m_map)
