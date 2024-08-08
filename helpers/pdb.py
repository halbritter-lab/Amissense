import numpy as np
import pandas as pd
from Bio import SeqIO, PDB
from Bio.PDB import PDBParser, PDBIO, DSSP
from pathlib import Path
from typing import Tuple


def extract_chain_id(uniprot_id: str, pdb_path: Path) -> str:
    """Fetches the chain ID for a given UniProt ID from the PDB file"""
    for record in SeqIO.parse(pdb_path, "pdb-seqres"):
        if uniprot_id in record.dbxrefs[0]:
            return record.annotations["chain"]
    raise KeyError(f"No matching chains found for {uniprot_id} in PDB file!")


def extract_secondary_structures(chain_id: str, pdb_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Extract secondary structure information from PDB file."""
    parser = PDBParser()
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]

    dssp = DSSP(model, str(pdb_path))
    max_residue = max(key[1][1] for key in dssp.property_dict if key[0] == chain_id)
    helices = np.zeros(max_residue + 1, dtype=int)
    sheets = np.zeros(max_residue + 1, dtype=int)

    for key, value in dssp.property_dict.items():
        if key[0] != chain_id:
            continue
        residue = key[1][1]
        ss = value[2]
        if ss == "H":  # Alpha-helix
            helices[residue] = 1
        elif ss == "E":  # Beta-sheet
            sheets[residue] = 1

    return helices, sheets


def extract_positional_confidences(chain_id: str, pdb_path: Path) -> np.ndarray:
    """Extract positional confidence information from PDB file."""
    parser = PDB.PDBParser()
    structure = parser.get_structure("protein", pdb_path)

    try:
        chain = structure[0][chain_id]
    except KeyError:
        raise ValueError(f"No chain found with ID: {chain_id}")

    max_residue = max(residue.id[1] for residue in chain if PDB.is_aa(residue))
    positional_confidences = np.zeros(shape=(max_residue + 1, 1))

    for residue in chain:
        if PDB.is_aa(residue) and "CA" in residue:
            sequence_number = residue.id[1]
            temperature_factor = residue["CA"].bfactor / 100.0
            positional_confidences[sequence_number] = temperature_factor

    return positional_confidences


def generate_pathogenicity_pdb(uniprot_id: str, predictions: pd.DataFrame, pdb_path: Path, out_dir: Path):
    """Generate a PDB file with pathogenicity values."""
    # Extract chain ID
    chain_id = extract_chain_id(uniprot_id, pdb_path)

    # Calculate average pathogenicity for each position
    predictions_grouped_means = predictions.groupby("protein_variant_pos")["pathogenicity"].mean()
    positional_means = predictions_grouped_means.reindex(
        range(0, predictions["protein_variant_pos"].max() + 1), fill_value=0
    ).to_numpy()

    print("Generating PDB with pathogenicity values")
    parser = PDBParser()
    structure = parser.get_structure("protein", str(pdb_path))
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.id[0] == " ":  # Check if it's a standard amino acid
                    residue.bfactor = positional_means[residue.id[1]]

    # Save new PDB file
    output_pdb_path = out_dir / f"{uniprot_id.upper()}_pathogenicity.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb_path))

    print(f"Generated PDB stored as {output_pdb_path}")


if __name__ == "__main__":
    pdb_path = Path("test/6lid.pdb")
    chain_id = extract_chain_id("P82251", pdb_path)
    pos_conf = extract_positional_confidences(chain_id, pdb_path)

    print(pos_conf)
    print(len(pos_conf))
