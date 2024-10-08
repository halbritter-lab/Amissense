import argparse
import logging
from pathlib import Path
import sys
import json
from Amissense.scripts.pipeline import run_pipeline
from Amissense.scripts.utils import get_uniprot_id, download_pdb, download_and_extract_alphamissense_predictions, load_config
from Amissense.version import __version__ as VERSION
from Amissense.scripts.json_generator import stream_tsv_file  # Import the function from json_generator.py

def setup_logging(log_level=logging.INFO, log_file=None):
    """Set up logging configuration."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file,
        filemode='a'  # Append mode for logging to a file
    )
    if not log_file:  # If no log file, also log to stdout
        logging.getLogger().addHandler(logging.StreamHandler())

def main():
    parser = argparse.ArgumentParser(description="Amissense CLI: Process AlphaMissense data and generate visualizations.")
    
    # Adding global flags with short and long options
    parser.add_argument('-l', '--log-level', help="Set the logging level", default="INFO")
    parser.add_argument('-f', '--log-file', help="Set the log output file", default=None)
    parser.add_argument('--config-path', type=Path, help="Path to the config.json file", default=None)
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {VERSION}')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand for running the full pipeline
    parser_pipeline = subparsers.add_parser("pipeline", help="Run the full Amissense pipeline.")
    parser_pipeline.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser_pipeline.add_argument('-u', '--uniprot-id', type=str, help="The UniProt ID of the protein (optional if gene-id is provided).")
    parser_pipeline.add_argument('-o', '--output-dir', type=str, default="out", help="Base directory to store output files.")
    parser_pipeline.add_argument('-e', '--experimental-pdb', type=str, default="", help="Path to experimental PDB file.")
    parser_pipeline.add_argument('-s', '--source', type=str, choices=['api', 'local'], default='api', help="Source for fetching AlphaMissense predictions (default: api).")
    parser_pipeline.add_argument('-i', '--organism-id', type=int, default=9606, help="Organism ID (default: 9606 for Homo sapiens).")
    parser_pipeline.add_argument('--include-clinvar-variants-in-pdb', action='store_true', help="Include ClinVar variants in the PDB visualizations.")

    # Subcommand for utils
    parser_utils = subparsers.add_parser("utils", help="Utility commands like fetching predictions, downloading PDB files, and querying UniProt.")
    utils_subparsers = parser_utils.add_subparsers(dest="utils_command", help="Available utility functions")

    # Subcommand for download_and_extract_alphamissense_predictions
    parser_download_predictions = utils_subparsers.add_parser("download-predictions", help="Download and extract AlphaMissense predictions.")
    parser_download_predictions.add_argument('-t', '--tmp-dir', type=Path, default=Path("tmp"), help="Temporary directory to download and extract the predictions.")
    
    # Subcommand for downloading PDB files
    parser_download_pdb = utils_subparsers.add_parser("download-pdb", help="Download a PDB file using its PDB ID.")
    parser_download_pdb.add_argument('-p', '--pdb-id', type=str, required=True, help="The PDB ID of the protein structure.")
    parser_download_pdb.add_argument('-o', '--output-dir', type=Path, default=Path("out/pdb"), help="Directory to save the downloaded PDB file.")

    # Subcommand for fetching UniProt ID
    parser_uniprot_query = utils_subparsers.add_parser("uniprot-query", help="Query UniProt for a gene's UniProt ID.")
    parser_uniprot_query.add_argument('-n', '--gene-name', type=str, required=True, help="The gene name to search for.")
    parser_uniprot_query.add_argument('-i', '--organism-id', type=int, required=True, help="The organism ID (e.g., 9606 for Homo sapiens).")

    # Subcommand for generating JSON files from TSV
    parser_json_generator = subparsers.add_parser("generate-json", help="Generate JSON files from AlphaMissense TSV data.")
    parser_json_generator.add_argument('tsv_file', type=Path, help="Path to the AlphaMissense_aa_substitutions.tsv.gz file.")
    parser_json_generator.add_argument('output_dir', type=Path, help="Directory to store generated JSON files.")
    
    # Parse arguments
    args = parser.parse_args()

    # Setup logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    setup_logging(log_level=log_level, log_file=args.log_file)

    # Load config with error handling
    try:
        config = load_config(args.config_path)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logging.critical(f"Failed to load configuration: {e}")
        sys.exit(1)

    # Handle subcommands
    if args.command == "pipeline":
        # Retrieve UniProt ID if not provided
        uniprot_id = args.uniprot_id
        if not uniprot_id:
            logging.info(f"Retrieving UniProt ID for gene {args.gene_id} and organism ID {args.organism_id}")
            uniprot_id = get_uniprot_id(args.gene_id, args.organism_id)
            if not uniprot_id:
                logging.error(f"Failed to retrieve UniProt ID for gene {args.gene_id}. Exiting.")
                sys.exit(1)
        
        run_pipeline(
            uniprot_id, 
            args.gene_id, 
            Path(args.output_dir), 
            Path(args.experimental_pdb) if args.experimental_pdb else None, 
            args.source, 
            include_clinvar_variants_in_pdb=args.include_clinvar_variants_in_pdb
        )

    elif args.command == "utils":
        # If no utils subcommand is given, show the help for utils subcommands
        if args.utils_command is None:
            parser_utils.print_help()
            sys.exit(1)

        if args.utils_command == "download-predictions":
            download_and_extract_alphamissense_predictions(args.tmp_dir)
        
        elif args.utils_command == "download-pdb":
            download_pdb(args.pdb_id, args.output_dir)
        
        elif args.utils_command == "uniprot-query":
            uniprot_id = get_uniprot_id(args.gene_name, args.organism_id)
            if uniprot_id:
                logging.info(f"UniProt ID: {uniprot_id}")
            else:
                logging.error("Failed to fetch UniProt ID.")
    
    elif args.command == "generate-json":
        stream_tsv_file(args.tsv_file, args.output_dir)

if __name__ == "__main__":
    main()
