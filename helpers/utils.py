import requests
import argparse
from pathlib import Path

def get_uniprot_id(gene_name: str, organism_id: int) -> str:
    """
    Query the UniProt REST API to get the primary accession (UniProt ID) for a given gene name and organism ID.
    The function returns the UniProt ID for the reviewed (Swiss-Prot) entry if available.

    Parameters:
    gene_name (str): The gene name to search for (e.g., "NAA10").
    organism_id (int): The organism ID (e.g., 9606 for Homo sapiens).

    Returns:
    str: The primary accession (UniProt ID) of the reviewed (Swiss-Prot) entry, or None if no reviewed entry is found.
    """
    # Define the base URL for the UniProt REST API
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Define the query parameters
    query = f"(organism_id:{organism_id}) AND (gene:{gene_name})"
    fields = "accession,reviewed,id,protein_name,gene_names,organism_name,length"
    params = {
        "query": query,
        "fields": fields,
        "format": "json"
    }

    try:
        # Send the GET request to the UniProt API
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an error for bad responses (4xx and 5xx)
        
        # Parse the JSON response
        data = response.json()

        # Iterate through the results to find the reviewed (Swiss-Prot) entry
        for entry in data.get("results", []):
            if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
                return entry.get("primaryAccession")

    except requests.RequestException as e:
        print(f"Error querying UniProt API: {e}")

    # Return None if no reviewed (Swiss-Prot) entry is found
    return None

def download_pdb_file(pdb_id: str, output_dir: Path) -> Path:
    """
    Download a PDB file from the RCSB PDB website using the given PDB ID and save it to the specified output directory.

    Parameters:
    pdb_id (str): The PDB ID of the protein structure (e.g., "6LID").
    output_dir (Path): The directory where the downloaded PDB file should be saved.

    Returns:
    Path: The path to the downloaded PDB file, or None if the download failed.
    """

    # Ensure the output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Construct the URL for the PDB file
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    try:
        # Send the GET request to download the PDB file
        response = requests.get(pdb_url)
        response.raise_for_status()  # Raise an error for bad responses (4xx and 5xx)

        # Define the output file path
        pdb_file_path = output_dir / f"{pdb_id}.pdb"

        # Write the content to the file
        with open(pdb_file_path, "wb") as file:
            file.write(response.content)
        
        print(f"Downloaded PDB file: {pdb_file_path}")
        return pdb_file_path

    except requests.RequestException as e:
        print(f"Error downloading PDB file {pdb_id}: {e}")
        return None

if __name__ == "__main__":
    # Setup argparse to accept command-line arguments
    parser = argparse.ArgumentParser(description="Utility script for querying UniProt and downloading PDB files.")
    
    # Subparsers for different commands (uniProt query and PDB download)
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subparser for querying UniProt
    uniprot_parser = subparsers.add_parser("uniprot", help="Query UniProt for a gene's UniProt ID based on gene name and organism ID")
    uniprot_parser.add_argument("gene_name", type=str, help="The gene name to search for (e.g., 'NAA10').")
    uniprot_parser.add_argument("organism_id", type=int, help="The organism ID (e.g., 9606 for Homo sapiens).")

    # Subparser for downloading PDB files
    pdb_parser = subparsers.add_parser("pdb", help="Download a PDB file using its PDB ID")
    pdb_parser.add_argument("pdb_id", type=str, help="The PDB ID of the protein structure (e.g., '6LID').")
    pdb_parser.add_argument("output_dir", type=Path, help="Directory to save the downloaded PDB file.")

    args = parser.parse_args()

    # Handle the subcommands
    if args.command == "uniprot":
        # Call the UniProt query function
        uniprot_id = get_uniprot_id(args.gene_name, args.organism_id)
        if uniprot_id:
            print(f"UniProt ID for {args.gene_name} in organism {args.organism_id}: {uniprot_id}")
        else:
            print(f"No reviewed UniProt ID found for {args.gene_name} in organism {args.organism_id}.")
    
    elif args.command == "pdb":
        # Call the PDB download function
        pdb_file_path = download_pdb_file(args.pdb_id, args.output_dir)
        if pdb_file_path:
            print(f"PDB file saved to: {pdb_file_path}")
        else:
            print("Failed to download the PDB file.")
