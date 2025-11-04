"""
ref boltz/src/boltz/data/parse/fasta.py
"""


import json, sys
from pathlib import Path

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

    > CHAIN_ID|ENTITY_TYPE|MSA_DIR
    SEQUENCE
    > CHAIN_ID|ENTITY_TYPE|MSA_DIR
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
            msg = f"Invalid record id: {seq_record.id}"
            raise ValueError(msg)

        header = seq_record.id.split("|")
        assert len(header) >= 2, f"Invalid record id: {seq_record.id}"

        #chain_id, entity_type = header[:2]
        entity_type,_ = header[:2]
        if entity_type.lower() not in {"protein", "dna", "rna", "ligand", "ion"}:
            msg = f"Invalid entity type: {entity_type}"
            raise ValueError(msg)
        #if chain_id == "":
        #    msg = "Empty chain id in input fasta!"
        #    raise ValueError(msg)
        if entity_type == "":
            msg = "Empty entity type in input fasta!"
            raise ValueError(msg)

    # Convert to json format
    sequences = []
    
    msa_ndx = 1
    for seq_record in records:
        # Get chain id, entity type and sequence
        header = seq_record.id.split("|")
        #chain_id, entity_type = header[:2]
        entity_type, _ = header[:2]
        entity_type = entity_type.upper()
        seq = str(seq_record.seq)

        if entity_type == "PROTEIN":
            msa_save_dir = Path(output_dir) / 'msa' / str(msa_ndx)
            if (len(header) == 3 and header[2] != "") or msa_save_dir.exists():
                assert (
                    entity_type.lower() == "protein"
                ), "MSA_DIR is only allowed for proteins"
                precomputed_msa_dir = str(msa_save_dir) if msa_save_dir.exists() else header[2]
                msa_id = {
                        "precomputed_msa_dir": precomputed_msa_dir,
                        "pairing_db": "uniref100"
                    }
            else:
                msa_id = {
                        "search_tool": "jackhmmer",
                        "pairing_db": "uniprot",
                        "pairing_db_fpath": UNIPROT_FA,
                        "non_pairing_db_fpath": MGNIFY_FA,
                        "msa_save_dir":  str(Path(output_dir) / 'msa')
                    }
                     #"msa_save_dir":  str(Path(output_dir) / 'msa' / str(msa_ndx) )
            msa_ndx += 1
            
            seq, ptm_datas = replace_and_track(seq)
            modifications = [{"ptmType": "CCD_" + ptm[0], "ptmPosition": ptm[1]} for  ptm in ptm_datas]
            
            molecule = {
                "proteinChain": {
                    "sequence": seq,
                    "count": 1,
                    "modifications": modifications,
                    "msa": msa_id,
                },
            }
        elif entity_type == "RNA":
            molecule = {
                "rnaSequence": {
                    "sequence": seq,
                    "count": 1
                },
            }
        elif entity_type == "DNA":
            molecule = {
                "dnaSequence": {
                    "sequence": seq,
                    "count": 1
                }
            }
        elif entity_type.upper() == "LIGAND":
            molecule = {
                "ligand": {
                    "ligand": seq,
                    "count": 1
                    
                }
            }
        elif entity_type.upper() == "ION":
            molecule = {
                "ion": {
                    "ion": seq,
                    "count": 1
                }
            }

        sequences.append(molecule)

    data = [{
        "name": "out",
        "sequences": sequences
    }]

    path = Path(output_dir)
    path.mkdir(exist_ok=True)
    jsonf = path / "inp.json"
    with open(jsonf, 'w') as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    parse_fasta(sys.argv[1], sys.argv[2])
