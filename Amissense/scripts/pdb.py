import numpy as np
import pandas as pd
from Bio import SeqIO, PDB
from Bio.PDB import PDBParser, PDBIO, DSSP
from pathlib import Path
from typing import Tuple
import argparse

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


def main():
    parser = argparse.ArgumentParser(description="Process PDB files and generate pathogenicity-modified PDB files.")
    
    # Add subcommands for different operations
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Subcommand for extracting chain ID
    parser_chain = subparsers.add_parser("extract_chain", help="Extract chain ID from PDB file based on UniProt ID.")
    parser_chain.add_argument("-u", "--uniprot-id", type=str, required=True, help="UniProt ID of the protein.")
    parser_chain.add_argument("-p", "--pdb-file", type=Path, required=True, help="Path to the PDB file.")
    
    # Subcommand for extracting secondary structures
    parser_secondary = subparsers.add_parser("extract_secondary", help="Extract secondary structures from PDB file.")
    parser_secondary.add_argument("-c", "--chain-id", type=str, required=True, help="Chain ID to extract structures from.")
    parser_secondary.add_argument("-p", "--pdb-file", type=Path, required=True, help="Path to the PDB file.")
    
    # Subcommand for extracting positional confidences
    parser_confidence = subparsers.add_parser("extract_confidence", help="Extract positional confidence from PDB file.")
    parser_confidence.add_argument("-c", "--chain-id", type=str, required=True, help="Chain ID to extract confidences from.")
    parser_confidence.add_argument("-p", "--pdb-file", type=Path, required=True, help="Path to the PDB file.")
    
    # Subcommand for generating pathogenicity PDB
    parser_generate = subparsers.add_parser("generate_pdb", help="Generate PDB with pathogenicity values.")
    parser_generate.add_argument("-u", "--uniprot-id", type=str, required=True, help="UniProt ID of the protein.")
    parser_generate.add_argument("-p", "--pdb-file", type=Path, required=True, help="Path to the PDB file.")
    parser_generate.add_argument("-o", "--out-dir", type=Path, required=True, help="Output directory for the new PDB file.")
    parser_generate.add_argument("-d", "--predictions-file", type=Path, required=True, help="Path to the predictions CSV file.")
    
    args = parser.parse_args()

    # Handle subcommands
    if args.command == "extract_chain":
        chain_id = extract_chain_id(args.uniprot_id, args.pdb_file)
        print(f"Extracted chain ID: {chain_id}")
    elif args.command == "extract_secondary":
        helices, sheets = extract_secondary_structures(args.chain_id, args.pdb_file)
        print("Helices:", helices)
        print("Sheets:", sheets)
    elif args.command == "extract_confidence":
        confidences = extract_positional_confidences(args.chain_id, args.pdb_file)
        print("Positional Confidences:", confidences)
    elif args.command == "generate_pdb":
        predictions = pd.read_csv(args.predictions_file)
        generate_pathogenicity_pdb(args.uniprot_id, predictions, args.pdb_file, args.out_dir)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
