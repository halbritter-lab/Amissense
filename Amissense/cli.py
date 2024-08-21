import argparse
import logging
from pathlib import Path
import sys
from Amissense.scripts.pipeline import run_pipeline  # Correct import
from Amissense.scripts.clinvar import fetch_clinvar_data
from Amissense.scripts.utils import get_uniprot_id, download_pdb_file
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
    
    # Subcommand for fetching ClinVar data
    parser_clinvar = subparsers.add_parser("clinvar", help="Fetch ClinVar data for a gene.")
    parser_clinvar.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID to fetch data for.")
    
    # Subcommand for downloading PDB files
    parser_pdb = subparsers.add_parser("pdb", help="Download a PDB file using its PDB ID.")
    parser_pdb.add_argument('-p', '--pdb-id', type=str, required=True, help="The PDB ID of the protein structure.")
    parser_pdb.add_argument('-o', '--output-dir', type=Path, required=True, help="Directory to save the downloaded PDB file.")
    
    # Subcommand for fetching UniProt ID
    parser_uniprot = subparsers.add_parser("uniprot", help="Query UniProt for a gene's UniProt ID.")
    parser_uniprot.add_argument('-n', '--gene-name', type=str, required=True, help="The gene name to search for.")
    parser_uniprot.add_argument('-i', '--organism-id', type=int, required=True, help="The organism ID (e.g., 9606 for Homo sapiens).")

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
    
    elif args.command == "clinvar":
        data = fetch_clinvar_data(args.gene_id)
        if data is not None:
            print(data)
        else:
            logging.error("Failed to fetch ClinVar data.")
    
    elif args.command == "pdb":
        download_pdb_file(args.pdb_id, args.output_dir)
    
    elif args.command == "uniprot":
        uniprot_id = get_uniprot_id(args.gene_name, args.organism_id)
        if uniprot_id:
            print(f"UniProt ID: {uniprot_id}")
        else:
            logging.error("Failed to fetch UniProt ID.")

if __name__ == "__main__":
    main()
