"""
ref boltz/src/boltz/data/parse/fasta.py
"""
import json, sys
from pathlib import Path
import string

from Bio import SeqIO

UNIPROT_FA = "/sugon_store/pub_data/databases/alphafold3/uniprot_all_2021_04.fa"
MGNIFY_FA  = "/sugon_store/pub_data/databases/alphafold3/mgy_clusters_2022_05.fa"

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


def parse_fasta(fasta_file, output_dir):
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

    Returns
    -------
    output_dir: Directory to save output.

    """
    # Read fasta file
    with Path(fasta_file).open("r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    # Make sure all records have a chain id and entity
    for seq_record in records:
        if "|" not in seq_record.id:
            count, entity_type = 1, 'protein'
            continue
            msg = f"Invalid record id: {seq_record.id}"
            raise ValueError(msg)

        header = seq_record.id.split("|")
        #assert len(header) >= 2, f"Invalid record id: {seq_record.id}"
        if header[0].isdigit():
            count, entity_type = header[:2]
        else:
            count, entity_type = 1, header[0]
        if entity_type.lower() not in {"protein", "dna", "rna", "ligand", "ion"}:
            msg = f"Invalid entity type: {entity_type}"
            raise ValueError(msg)
        if entity_type == "":
            msg = "Empty entity type in input fasta!"
            raise ValueError(msg)
    
    # Convert to json format
    sequences = []
    
    #msa_ndx = 1
    for ndx, seq_record in enumerate(records):
        if "|" not in seq_record.id:
            count, entity_type = 1, 'protein'
        else:
            # Get chain id, entity type and sequence
            header = seq_record.id.split("|")
            if header[0].isdigit():
                count, entity_type = header[:2]
            else:
                count, entity_type = 1, header[0]
        count = int(count)
        entity_type = entity_type.upper()
        seq = str(seq_record.seq)

        if entity_type == "PROTEIN":
            
            seq, ptm_datas = replace_and_track(seq)
            modifications = [{"ptmType": ptm[0], "ptmPosition": ptm[1]} for  ptm in ptm_datas]
            
            molecule = {
                "protein": {
                "id": string.ascii_uppercase[ndx],
                "sequence": seq,
                "modifications": modifications
                }
            }
        elif entity_type == "RNA":
             molecule = {
                "rna": {
                "id": string.ascii_uppercase[ndx],
                "sequence": seq,
                "modifications": []
                }
            }
        elif entity_type == "DNA":
            molecule = {
                "dna": {
                "id": string.ascii_uppercase[ndx],
                "sequence": seq,
                "modifications": []
                }
            }
        elif entity_type.upper() == "LIGAND":
            
            molecule = {
                "ligand": {
                "id": string.ascii_uppercase[ndx],
                }
            }
            #if seq.startswith('['):
            #    molecule['ligand']['ccdCodes'] = eval(seq)
            if len(seq) == 3:
                molecule['ligand']['ccdCodes'] = [seq]
            else:
                molecule['ligand']['smiles'] = seq
                
        elif entity_type.upper() == "ION":
            molecule = {
                "ion": {
                    "ion": seq,
                    "count": count
                }
            }

        sequences.append(molecule)
    
    data = {
        "name": "out",
        "modelSeeds": [1],  # At least one seed required.
        "sequences": sequences,
        "dialect": "alphafold3",  # Required
        "version": 2  # Required
        }

    path = Path(output_dir)
    path.mkdir(exist_ok=True)
    jsonf = path / "inp.json"
    with open(jsonf, 'w') as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    parse_fasta(sys.argv[1], sys.argv[2])
