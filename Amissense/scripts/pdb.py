import numpy as np
import pandas as pd
from Bio import SeqIO, PDB
from Bio.PDB import PDBParser, PDBIO, DSSP
from pathlib import Path
from typing import Tuple
import argparse
import logging
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning  # Import the specific warning
import Amissense.scripts.utils as utils

# Load configuration from utils module
try:
    config = utils.load_config()
except SystemExit:
    logging.critical("Unable to load configuration. Exiting.")
    raise

# Setup logging configuration
logging.basicConfig(
    level=config["logging"]["default_level"],
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)

# Suppress specific PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

# Load Chimera color ranges from config.json
chimera_color_ranges = config.get("chimera_color_ranges", {})

def extract_chain_id(uniprot_id: str, pdb_path: Path) -> str:
    """Fetches the chain ID for a given UniProt ID from the PDB file"""
    try:
        for record in SeqIO.parse(pdb_path, "pdb-seqres"):
            if uniprot_id in record.dbxrefs[0]:
                return record.annotations["chain"]
        logging.error(f"No matching chains found for {uniprot_id} in PDB file!")
        raise KeyError(f"No matching chains found for {uniprot_id} in PDB file!")
    except FileNotFoundError:
        logging.error(f"PDB file not found: {pdb_path}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error in extracting chain ID: {e}")
        raise

def extract_secondary_structures(chain_id: str, pdb_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Extract secondary structure information from PDB file."""
    try:
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
    except FileNotFoundError:
        logging.error(f"PDB file not found: {pdb_path}")
        raise
    except KeyError:
        logging.error(f"Chain ID {chain_id} not found in PDB file: {pdb_path}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error in extracting secondary structures: {e}")
        raise

def extract_positional_confidences(chain_id: str, pdb_path: Path) -> np.ndarray:
    """Extract positional confidence information from PDB file."""
    try:
        parser = PDB.PDBParser()
        structure = parser.get_structure("protein", pdb_path)

        try:
            chain = structure[0][chain_id]
        except KeyError:
            logging.error(f"No chain found with ID: {chain_id}")
            raise ValueError(f"No chain found with ID: {chain_id}")

        max_residue = max(residue.id[1] for residue in chain if PDB.is_aa(residue))
        positional_confidences = np.zeros(shape=(max_residue + 1, 1))

        for residue in chain:
            if PDB.is_aa(residue) and "CA" in residue:
                sequence_number = residue.id[1]
                temperature_factor = residue["CA"].bfactor / 100.0
                positional_confidences[sequence_number] = temperature_factor

        return positional_confidences
    except FileNotFoundError:
        logging.error(f"PDB file not found: {pdb_path}")
        raise
    except KeyError as e:
        logging.error(f"Key error in extracting positional confidences: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error in extracting positional confidences: {e}")
        raise

def generate_pymol_script(uniprot_id: str, predictions_grouped_means: pd.Series, pdb_path: Path, out_dir: Path, clinvar_variants: pd.Series = None):
    """
    Generate a PyMOL script that colors residues based on pathogenicity values and shows ClinVar variants as spheres.
    """
    try:
        # Create the output directory if it doesn't exist
        out_dir.mkdir(parents=True, exist_ok=True)

        # Path to save the PyMOL script
        script_path = out_dir / f"{uniprot_id}_pymol_script.pml"

        # Determine if we need to use a relative path for the PDB file
        pdb_relative_path = pdb_path.relative_to(out_dir) if pdb_path.parent == out_dir else pdb_path

        # Open the script file for writing
        with open(script_path, "w") as script_file:
            # PyMOL script header
            script_file.write("# PyMOL script to color residues based on pathogenicity values and show ClinVar variants as spheres\n")
            script_file.write("from pymol import cmd\n")
            
            # Load the PDB file
            script_file.write(f"cmd.load('{pdb_relative_path}')\n")
            
            # Iterate over predictions and generate PyMOL commands for coloring
            for residue_pos, pathogenicity in predictions_grouped_means.items():
                color = get_color_from_config(pathogenicity)
                script_file.write(f"cmd.color('{color}', 'resi {residue_pos}')\n")

            # Optionally show ClinVar variants as spheres
            if clinvar_variants is not None:
                for residue_pos, score in clinvar_variants.items():
                    script_file.write(f"cmd.show('spheres', 'resi {residue_pos}')\n")

            # Save the session
            session_filename = f"{uniprot_id}_pathogenicity_session.pse"
            script_file.write(f"cmd.save('{session_filename}')\n")

        logging.info(f"PyMOL script saved at {script_path}")

    except Exception as e:
        logging.error(f"Error generating PyMOL script: {e}")
        raise

def generate_chimera_session(uniprot_id: str, predictions_grouped_means: pd.Series, pdb_path: Path, out_dir: Path, clinvar_variants: pd.Series = None):
    """
    Generate a Chimera session script that colors residues based on pathogenicity values and shows ClinVar variants as spheres.

    Parameters:
    - uniprot_id: The UniProt ID of the protein.
    - predictions_grouped_means: A pandas Series where the index is the residue number, and the value is the pathogenicity score.
    - pdb_path: The path to the PDB file to load in Chimera.
    - out_dir: The directory where the Chimera session script will be saved.
    - clinvar_variants: Optional pandas Series where the index is the residue number, and the value indicates the pathogenicity score for ClinVar variants.
    """
    try:
        # Create the output directory if it doesn't exist
        out_dir.mkdir(parents=True, exist_ok=True)

        # Path to save the Chimera session script
        script_path = out_dir / f"{uniprot_id}_chimera_session.py"

        # Determine if we need to use a relative path for the PDB file
        pdb_relative_path = pdb_path.relative_to(out_dir) if pdb_path.parent == out_dir else pdb_path

        # Open the script file for writing
        with open(script_path, "w") as script_file:
            # Chimera script header
            script_file.write("# Chimera script to color residues based on pathogenicity values and show ClinVar variants as spheres\n")
            script_file.write("import chimera\n")
            script_file.write("from chimera import runCommand\n")
            
            # Load the PDB file
            script_file.write(f"runCommand('open {pdb_relative_path}')\n")
            
            # Iterate over predictions and generate Chimera commands for coloring
            for residue_pos, pathogenicity in predictions_grouped_means.items():
                # Determine the color based on the score ranges defined in the config
                color = get_color_from_config(pathogenicity)
                script_file.write(f"runCommand('color {color} :{residue_pos}')\n")

            # Optionally show ClinVar variants as spheres
            if clinvar_variants is not None:
                for residue_pos, score in clinvar_variants.items():
                    # Ensure the residue is displayed and then represented as a sphere
                    script_file.write(f"runCommand('display :{residue_pos}')\n")
                    script_file.write(f"runCommand('represent sphere :{residue_pos}')\n")

            # Save the session (use only the filename, not the full path)
            session_filename = f"{uniprot_id}_pathogenicity_session.py"
            script_file.write(f"runCommand('save {session_filename}')\n")

        logging.info(f"Chimera session script saved at {script_path}")

    except Exception as e:
        logging.error(f"Error generating Chimera session: {e}")
        raise

def get_color_from_config(score):
    """
    Map a pathogenicity score to a color using the chimera_color_ranges in config.json.
    """
    for range_key, color in chimera_color_ranges.items():
        low, high = map(float, range_key.split("-"))
        if low <= score <= high:
            return color
    logging.warning(f"No color found for score {score}. Please check the configured ranges in config.json. Using default 'gray'.")
    return "gray"

def generate_pathogenicity_pdb(uniprot_id: str, predictions: pd.DataFrame, pdb_path: Path, out_dir: Path):
    """Generate a PDB file with pathogenicity values."""
    try:
        # Extract chain ID
        chain_id = extract_chain_id(uniprot_id, pdb_path)

        # Calculate average pathogenicity for each position
        predictions_grouped_means = predictions.groupby("protein_variant_pos")["pathogenicity"].mean()
        positional_means = predictions_grouped_means.reindex(
            range(0, predictions["protein_variant_pos"].max() + 1), fill_value=0
        ).to_numpy()

        logging.info("Generating PDB with pathogenicity values")
        parser = PDBParser()
        structure = parser.get_structure("protein", str(pdb_path))
        atom_list = structure.get_atoms()
        for atom in atom_list:
            if atom.full_id[2] != chain_id:
                continue
            if atom.full_id[3][0] != " ":
                continue
            atom.set_bfactor(positional_means[atom.full_id[3][1]])

        # Save new PDB file
        output_pdb_path = out_dir / f"{uniprot_id.upper()}_pathogenicity.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb_path))

        logging.info(f"Generated PDB stored as {output_pdb_path}")

    except FileNotFoundError:
        logging.error(f"PDB file not found: {pdb_path}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error in generating pathogenicity PDB: {e}")
        raise

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
    try:
        if args.command == "extract_chain":
            chain_id = extract_chain_id(args.uniprot_id, args.pdb_file)
            logging.info(f"Extracted chain ID: {chain_id}")
        elif args.command == "extract_secondary":
            helices, sheets = extract_secondary_structures(args.chain_id, args.pdb_file)
            logging.info(f"Helices: {helices}")
            logging.info(f"Sheets: {sheets}")
        elif args.command == "extract_confidence":
            confidences = extract_positional_confidences(args.chain_id, args.pdb_file)
            logging.info(f"Positional Confidences: {confidences}")
        elif args.command == "generate_pdb":
            predictions = pd.read_csv(args.predictions_file)
            generate_pathogenicity_pdb(args.uniprot_id, predictions, args.pdb_file, args.out_dir)
        else:
            parser.print_help()
    except Exception as e:
        logging.critical(f"Failed to process command: {e}")
        raise

if __name__ == "__main__":
    main()
