import argparse
import requests
import pandas as pd
from pathlib import Path

import Amissense.scripts.pdb as pdb_module
import Amissense.scripts.graphs as graphs_module
import Amissense.scripts.clinvar as clinvar_module
import Amissense.scripts.utils as utils

TMP_DIR = Path("tmp")

def ensure_directory_exist(dir: Path):
    """Create directories if they don't exist."""
    dir.mkdir(parents=True, exist_ok=True)

def download_alphafold_pdb(uniprot_id: str) -> Path:
    """Fetch AlphaFold PDB file for a given UniProt ID."""
    download_dir = TMP_DIR / "pdb_files"
    ensure_directory_exist(download_dir)
    alphafold_pdb_path = download_dir / f"{uniprot_id.upper()}_alphafold.pdb"

    if not alphafold_pdb_path.exists():
        alphafold_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"
        print(f"Downloading PDB file from {alphafold_api}")
        response = requests.get(alphafold_api, timeout=15)
        response.raise_for_status()
        data = response.json()
        pdb_url = data[0]["pdbUrl"]
        pdb_response = requests.get(pdb_url)
        pdb_response.raise_for_status()
        alphafold_pdb_path.write_bytes(pdb_response.content)
        print(f"Download completed and saved to {alphafold_pdb_path}")

    return alphafold_pdb_path

def run_pipeline(uniprot_id: str, gene_id: str, output_dir: Path, experimental_pdb: Path):
    """Main function to run the pipeline."""
    ensure_directory_exist(output_dir)

    try:
        # Fetch AlphaMissense predictions from JSON
        json_dir = TMP_DIR / "json"
        predictions = utils.get_predictions_from_json_for_uniprot_id(uniprot_id, json_dir)
        predictions.to_csv(output_dir / f"{uniprot_id}_AM_pathogenicity_predictions.csv", index=False)

        # Use experimental PDB if provided, otherwise download from AlphaFold
        if experimental_pdb == Path():
            pdb_path = download_alphafold_pdb(uniprot_id)
        else:
            print("Using experimental PDB...")
            pdb_path = experimental_pdb

        # Extract PDB details
        chain_id = pdb_module.extract_chain_id(uniprot_id, pdb_path)
        helices, sheets = pdb_module.extract_secondary_structures(chain_id, pdb_path)
        pdb_confidences = (
            pdb_module.extract_positional_confidences(chain_id, pdb_path) if experimental_pdb == Path() else None
        )

        # Generate PDB with pathogenicity and create visualizations
        pdb_module.generate_pathogenicity_pdb(uniprot_id, predictions, pdb_path, output_dir)
        graphs_module.plot_predictions_heatmap(uniprot_id, predictions, output_dir)
        graphs_module.plot_predictions_line_graph(uniprot_id, predictions, output_dir, helices, sheets, pdb_confidences)

        # Fetch and merge ClinVar data, then create visualizations
        clinvar_data = clinvar_module.fetch_clinvar_data(gene_id)
        clinvar_merged_data = clinvar_module.merge_missense_data(clinvar_data, predictions)
        clinvar_merged_data.to_csv(str(output_dir / f"{gene_id}_clinvar_AM.csv"), index=False)

        graphs_module.plot_clinvar_scatter(gene_id, predictions, clinvar_merged_data, output_dir)
        graphs_module.plot_clinvar_sankey(gene_id, predictions, clinvar_merged_data, output_dir)

    except requests.HTTPError as http_err:
        print(f"Unexpected error during download: {http_err}")
    except KeyError as key_err:
        print(key_err)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the AlphaMissense data processing pipeline.")
    parser.add_argument('-u', '--uniprot-id', type=str, required=True, help="The UniProt ID of the protein.")
    parser.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser.add_argument('-o', '--output-dir', type=str, default="out", help="Directory to store output files.")
    parser.add_argument('-e', '--experimental-pdb', type=str, default="", help="Path to experimental PDB file.")
    
    args = parser.parse_args()
    run_pipeline(args.uniprot_id, args.gene_id, Path(args.output_dir), Path(args.experimental_pdb))
