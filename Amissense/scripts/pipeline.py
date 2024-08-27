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
from Amissense.scripts.pdb import generate_chimera_session

# Load configuration with error handling
try:
    config = utils.load_config()
except SystemExit:
    logging.critical("Unable to load configuration. Exiting.")
    raise

def run_pipeline(uniprot_id: str, gene_id: str, base_output_dir: Path, experimental_pdb: Path = None, source: str = 'api'):
    """Main function to run the pipeline."""
    
    # Ensure the base directory exists
    try:
        output_dir = utils.generate_output_directory(base_output_dir, gene_id, uniprot_id)

        # Create subdirectories
        figures_dir = output_dir / "figures"
        tables_dir = output_dir / "tables"
        pdb_dir = output_dir / "pdb"
        
        utils.ensure_directory_exists(figures_dir)
        utils.ensure_directory_exists(tables_dir)
        utils.ensure_directory_exists(pdb_dir)
    except OSError as e:
        logging.critical(f"Failed to create output directories: {e}")
        return

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

        # Check for experimental PDB or download from AlphaFold/RCSB
        if experimental_pdb and not experimental_pdb.exists():
            # Extract PDB ID from the provided experimental PDB file path name (assumes filename contains PDB ID)
            pdb_id = experimental_pdb.stem[:4]  # Adjust based on actual naming convention
            logging.warning(f"Experimental PDB file not found. Attempting to download PDB ID: {pdb_id}")
            pdb_path = utils.download_pdb(uniprot_id, pdb_dir, pdb_id=pdb_id)
        elif experimental_pdb and experimental_pdb.exists():
            pdb_path = experimental_pdb
        else:
            pdb_path = utils.download_pdb(uniprot_id, pdb_dir)

        if not pdb_path:
            logging.error("Failed to obtain PDB file. Pipeline cannot proceed.")
            return

        # Copy the original PDB file to the output directory
        original_pdb_path = pdb_dir / f"{uniprot_id.upper()}_original.pdb"
        pdb_path.rename(original_pdb_path)

        # Extract PDB details
        chain_id = pdb_module.extract_chain_id(uniprot_id, original_pdb_path)
        helices, sheets = pdb_module.extract_secondary_structures(chain_id, original_pdb_path)
        pdb_confidences = pdb_module.extract_positional_confidences(chain_id, original_pdb_path) if not experimental_pdb else None

        # Generate PDB with pathogenicity and create visualizations
        pdb_module.generate_pathogenicity_pdb(uniprot_id, predictions, original_pdb_path, pdb_dir)
        graphs_module.plot_predictions_heatmap(uniprot_id, predictions, figures_dir)
        graphs_module.plot_predictions_line_graph(uniprot_id, predictions, figures_dir, helices, sheets, pdb_confidences)

        # Generate Chimera session
        predictions_grouped_means = predictions.groupby("protein_variant_pos")["pathogenicity"].mean()
        generate_chimera_session(uniprot_id, predictions_grouped_means, original_pdb_path, pdb_dir)

        # Fetch and merge ClinVar data, then create visualizations
        clinvar_data = clinvar_module.fetch_clinvar_data(gene_id)
        clinvar_merged_data = clinvar_module.merge_missense_data(clinvar_data, predictions)
        
        clinvar_file = tables_dir / f"{gene_id}_{uniprot_id}_clinvar_AM_{date_str}.csv"
        clinvar_merged_data.to_csv(clinvar_file, index=False)

        graphs_module.plot_clinvar_scatter(gene_id, predictions, clinvar_merged_data, figures_dir)
        graphs_module.plot_clinvar_sankey(gene_id, predictions, clinvar_merged_data, figures_dir)

    except requests.RequestException as e:
        logging.error(f"Network-related error: {e}")
    except FileNotFoundError as e:
        logging.error(f"File-related error: {e}")
    except KeyError as e:
        logging.error(f"Key-related error: {e}")
    except pd.errors.EmptyDataError as e:
        logging.error(f"Data-related error (empty file): {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the AlphaMissense data processing pipeline.")
    parser.add_argument('-u', '--uniprot-id', type=str, required=True, help="The UniProt ID of the protein.")
    parser.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser.add_argument('-o', '--output-dir', type=str, default=config['directories']['output_dir'], help="Base directory to store output files.")
    parser.add_argument('-e', '--experimental-pdb', type=str, default=None, help="Path to experimental PDB file. If provided, this will be used instead of downloading.")
    parser.add_argument('-s', '--source', type=str, choices=['api', 'local'], default='api', help="Source for fetching AlphaMissense predictions (default: api).")

    args = parser.parse_args()

    try:
        run_pipeline(args.uniprot_id, args.gene_id, Path(args.output_dir), Path(args.experimental_pdb) if args.experimental_pdb else None, args.source)
    except Exception as e:
        logging.critical(f"Pipeline failed: {e}")
        raise
