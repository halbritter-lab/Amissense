import argparse
import logging
from pathlib import Path
import sys
from Amissense.scripts.pipeline import run_pipeline  # Correct import
from Amissense.scripts.utils import get_uniprot_id, download_pdb_file, download_and_extract_alphamissense_predictions
from Amissense.version import __version__ as VERSION  # Assuming version is stored in Amissense/version.py

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
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {VERSION}')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand for running the full pipeline
    parser_pipeline = subparsers.add_parser("pipeline", help="Run the full Amissense pipeline.")
    parser_pipeline.add_argument('-u', '--uniprot-id', type=str, required=True, help="The UniProt ID of the protein.")
    parser_pipeline.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser_pipeline.add_argument('-o', '--output-dir', type=str, default="out", help="Directory to store output files.")
    parser_pipeline.add_argument('-e', '--experimental-pdb', type=str, default="", help="Path to experimental PDB file.")
    
    # Subcommand for utils
    parser_utils = subparsers.add_parser("utils", help="Utility commands like fetching predictions, downloading PDB files, and querying UniProt.")
    utils_subparsers = parser_utils.add_subparsers(dest="utils_command", help="Available utility functions")

    # Subcommand for download_and_extract_alphamissense_predictions
    parser_download_predictions = utils_subparsers.add_parser("download-predictions", help="Download and extract AlphaMissense predictions.")
    parser_download_predictions.add_argument('-t', '--tmp-dir', type=Path, required=True, help="Temporary directory to download and extract the predictions.")
    
    # Subcommand for downloading PDB files
    parser_download_pdb = utils_subparsers.add_parser("download-pdb", help="Download a PDB file using its PDB ID.")
    parser_download_pdb.add_argument('-p', '--pdb-id', type=str, required=True, help="The PDB ID of the protein structure.")
    parser_download_pdb.add_argument('-o', '--output-dir', type=Path, required=True, help="Directory to save the downloaded PDB file.")

    # Subcommand for fetching UniProt ID
    parser_uniprot_query = utils_subparsers.add_parser("uniprot-query", help="Query UniProt for a gene's UniProt ID.")
    parser_uniprot_query.add_argument('-n', '--gene-name', type=str, required=True, help="The gene name to search for.")
    parser_uniprot_query.add_argument('-i', '--organism-id', type=int, required=True, help="The organism ID (e.g., 9606 for Homo sapiens).")

    args = parser.parse_args()

    # If no command is given, print the help message
    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Setup logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    setup_logging(log_level=log_level, log_file=args.log_file)

    # Handle subcommands
    if args.command == "pipeline":
        run_pipeline(args.uniprot_id, args.gene_id, Path(args.output_dir), Path(args.experimental_pdb))
    
    elif args.command == "utils":
        # If no utils subcommand is given, show the help for utils subcommands
        if args.utils_command is None:
            parser_utils.print_help()
            sys.exit(1)

        if args.utils_command == "download-predictions":
            download_and_extract_alphamissense_predictions(args.tmp_dir)
        
        elif args.utils_command == "download-pdb":
            download_pdb_file(args.pdb_id, args.output_dir)
        
        elif args.utils_command == "uniprot-query":
            uniprot_id = get_uniprot_id(args.gene_name, args.organism_id)
            if uniprot_id:
                print(f"UniProt ID: {uniprot_id}")
            else:
                logging.error("Failed to fetch UniProt ID.")

if __name__ == "__main__":
    main()
