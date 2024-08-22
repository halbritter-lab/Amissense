import argparse
import logging
import requests
import pandas as pd
from pathlib import Path

import Amissense.scripts.pdb as pdb_module
import Amissense.scripts.graphs as graphs_module
import Amissense.scripts.clinvar as clinvar_module
import Amissense.scripts.utils as utils

# Load configuration from config.json
config = utils.load_config()

# Use the TMP_DIR from the config file
TMP_DIR = Path(config['directories']['tmp_dir'])

def run_pipeline(uniprot_id: str, gene_id: str, output_dir: Path, experimental_pdb: Path, source: str):
    """Main function to run the pipeline."""
    utils.ensure_directory_exists(output_dir)

    try:
        # Fetch AlphaMissense predictions based on the source
        if source == 'api':
            predictions = utils.get_predictions_from_static_api(uniprot_id)
        else:
            json_dir = Path(config['directories']['json_dir'])
            predictions = utils.get_predictions_from_json_for_uniprot_id(uniprot_id, json_dir)
        
        predictions.to_csv(output_dir / f"{uniprot_id}_AM_pathogenicity_predictions.csv", index=False)

        # Use experimental PDB if provided, otherwise download from AlphaFold
        pdb_path = experimental_pdb if experimental_pdb else download_alphafold_pdb(uniprot_id)

        # Extract PDB details
        chain_id = pdb_module.extract_chain_id(uniprot_id, pdb_path)
        helices, sheets = pdb_module.extract_secondary_structures(chain_id, pdb_path)
        pdb_confidences = pdb_module.extract_positional_confidences(chain_id, pdb_path) if not experimental_pdb else None

        # Generate PDB with pathogenicity and create visualizations
        pdb_module.generate_pathogenicity_pdb(uniprot_id, predictions, pdb_path, output_dir)
        graphs_module.plot_predictions_heatmap(uniprot_id, predictions, output_dir)
        graphs_module.plot_predictions_line_graph(uniprot_id, predictions, output_dir, helices, sheets, pdb_confidences)

        # Fetch and merge ClinVar data, then create visualizations
        clinvar_data = clinvar_module.fetch_clinvar_data(gene_id)
        clinvar_merged_data = clinvar_module.merge_missense_data(clinvar_data, predictions)
        clinvar_merged_data.to_csv(output_dir / f"{gene_id}_clinvar_AM.csv", index=False)

        graphs_module.plot_clinvar_scatter(gene_id, predictions, clinvar_merged_data, output_dir)
        graphs_module.plot_clinvar_sankey(gene_id, predictions, clinvar_merged_data, output_dir)

    except requests.HTTPError as http_err:
        logging.error(f"Unexpected error during download: {http_err}")
    except KeyError as key_err:
        logging.error(f"Key error: {key_err}")

def download_alphafold_pdb(uniprot_id: str) -> Path:
    """Fetch AlphaFold PDB file for a given UniProt ID."""
    download_dir = Path(config['directories']['pdb_dir'])
    utils.ensure_directory_exists(download_dir)
    alphafold_pdb_path = download_dir / f"{uniprot_id.upper()}_alphafold.pdb"

    if not alphafold_pdb_path.exists():
        alphafold_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"
        logging.info(f"Downloading PDB file from {alphafold_api}")
        response = requests.get(alphafold_api, timeout=15)
        response.raise_for_status()
        data = response.json()
        pdb_url = data[0]["pdbUrl"]
        pdb_response = requests.get(pdb_url)
        pdb_response.raise_for_status()
        alphafold_pdb_path.write_bytes(pdb_response.content)
        logging.info(f"Download completed and saved to {alphafold_pdb_path}")

    return alphafold_pdb_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the AlphaMissense data processing pipeline.")
    parser.add_argument('-u', '--uniprot-id', type=str, required=True, help="The UniProt ID of the protein.")
    parser.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser.add_argument('-o', '--output-dir', type=str, default=config['directories']['output_dir'], help="Directory to store output files.")
    parser.add_argument('-e', '--experimental-pdb', type=str, default="", help="Path to experimental PDB file.")
    parser.add_argument('-s', '--source', type=str, choices=['api', 'local'], default='api', help="Source for fetching AlphaMissense predictions (default: api).")

    args = parser.parse_args()
    run_pipeline(args.uniprot_id, args.gene_id, Path(args.output_dir), Path(args.experimental_pdb), args.source)
