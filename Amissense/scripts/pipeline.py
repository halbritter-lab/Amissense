import argparse
import logging
import requests
import pandas as pd
from pathlib import Path
from datetime import datetime

import Amissense.scripts.pdb as pdb_module
import Amissense.scripts.graphs as graphs_module
import Amissense.scripts.clinvar as clinvar_module
import Amissense.scripts.utils as utils

# Load configuration from config.json
config = utils.load_config()

def run_pipeline(uniprot_id: str, gene_id: str, base_output_dir: Path, experimental_pdb: Path = None, source: str = 'api'):
    """Main function to run the pipeline."""
    # Generate structured output directory name with date
    output_dir = utils.generate_output_directory(base_output_dir, gene_id, uniprot_id)

    # Create subdirectories for figures, tables, and pdb files directly within the output directory
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    pdb_dir = output_dir / "pdb"
    
    utils.ensure_directory_exists(figures_dir)
    utils.ensure_directory_exists(tables_dir)
    utils.ensure_directory_exists(pdb_dir)

    date_str = datetime.now().strftime("%Y-%m-%d")

    try:
        # Fetch AlphaMissense predictions based on the source
        if source == 'api':
            predictions = utils.get_predictions_from_static_api(uniprot_id)
        else:
            json_dir = Path(config['directories']['json_dir'])
            predictions = utils.get_predictions_from_json_for_uniprot_id(uniprot_id, json_dir)
        
        # Save predictions to the tables directory with date
        predictions_file = tables_dir / f"{gene_id}_{uniprot_id}_AM_pathogenicity_predictions_{date_str}.csv"
        predictions.to_csv(predictions_file, index=False)

        # Use experimental PDB if provided, otherwise download from AlphaFold
        pdb_path = experimental_pdb if experimental_pdb else utils.download_pdb(uniprot_id, pdb_dir)

        # Extract PDB details
        chain_id = pdb_module.extract_chain_id(uniprot_id, pdb_path)
        helices, sheets = pdb_module.extract_secondary_structures(chain_id, pdb_path)
        pdb_confidences = pdb_module.extract_positional_confidences(chain_id, pdb_path) if not experimental_pdb else None

        # Generate PDB with pathogenicity and create visualizations
        pdb_module.generate_pathogenicity_pdb(uniprot_id, predictions, pdb_path, pdb_dir)
        graphs_module.plot_predictions_heatmap(uniprot_id, predictions, figures_dir)
        graphs_module.plot_predictions_line_graph(uniprot_id, predictions, figures_dir, helices, sheets, pdb_confidences)

        # Fetch and merge ClinVar data, then create visualizations
        clinvar_data = clinvar_module.fetch_clinvar_data(gene_id)
        clinvar_merged_data = clinvar_module.merge_missense_data(clinvar_data, predictions)
        
        clinvar_file = tables_dir / f"{gene_id}_{uniprot_id}_clinvar_AM_{date_str}.csv"
        clinvar_merged_data.to_csv(clinvar_file, index=False)

        graphs_module.plot_clinvar_scatter(gene_id, predictions, clinvar_merged_data, figures_dir)
        graphs_module.plot_clinvar_sankey(gene_id, predictions, clinvar_merged_data, figures_dir)

    except requests.HTTPError as http_err:
        logging.error(f"Unexpected error during download: {http_err}")
    except KeyError as key_err:
        logging.error(f"Key error: {key_err}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the AlphaMissense data processing pipeline.")
    parser.add_argument('-u', '--uniprot-id', type=str, required=True, help="The UniProt ID of the protein.")
    parser.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser.add_argument('-o', '--output-dir', type=str, default=config['directories']['output_dir'], help="Base directory to store output files.")
    parser.add_argument('-e', '--experimental-pdb', type=str, default=None, help="Path to experimental PDB file. If provided, this will be used instead of downloading.")
    parser.add_argument('-s', '--source', type=str, choices=['api', 'local'], default='api', help="Source for fetching AlphaMissense predictions (default: api).")

    args = parser.parse_args()
    run_pipeline(args.uniprot_id, args.gene_id, Path(args.output_dir), Path(args.experimental_pdb) if args.experimental_pdb else None, args.source)
