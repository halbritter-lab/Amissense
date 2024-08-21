import requests
import gzip
import shutil
import pandas as pd
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
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"(organism_id:{organism_id}) AND (gene:{gene_name})"
    fields = "accession,reviewed,id,protein_name,gene_names,organism_name,length"
    params = {"query": query, "fields": fields, "format": "json"}

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        for entry in data.get("results", []):
            if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
                return entry.get("primaryAccession")
    except requests.RequestException as e:
        print(f"Error querying UniProt API: {e}")

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
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    try:
        response = requests.get(pdb_url)
        response.raise_for_status()
        pdb_file_path = output_dir / f"{pdb_id}.pdb"
        with open(pdb_file_path, "wb") as file:
            file.write(response.content)
        print(f"Downloaded PDB file: {pdb_file_path}")
        return pdb_file_path
    except requests.RequestException as e:
        print(f"Error downloading PDB file {pdb_id}: {e}")
        return None

def download_and_extract_alphamissense_predictions(tmp_dir: Path) -> Path:
    """
    Download and extract AlphaMissense predictions.

    Parameters:
    tmp_dir (Path): The temporary directory where the prediction file will be downloaded and extracted.

    Returns:
    Path: The path to the extracted TSV file.
    """
    tmp_dir.mkdir(parents=True, exist_ok=True)
    url = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
    file_path = tmp_dir / Path(url).name
    tsv_path = file_path.with_suffix("")

    if not tsv_path.exists():
        if not file_path.exists():
            print(f"Downloading AlphaMissense predictions from {url}")
            response = requests.get(url)
            response.raise_for_status()
            file_path.write_bytes(response.content)
            print("Download completed!")
        
        print("Extracting predictions")
        with gzip.open(file_path, "rb") as f_in, tsv_path.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print("Extraction completed!")
    
    return tsv_path

def get_predictions_from_am_tsv_for_uniprot_id(uniprot_id: str, missense_tsv: Path) -> pd.DataFrame:
    """
    Fetch AlphaMissense predictions for a specific UniProt ID.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    missense_tsv (Path): The path to the AlphaMissense TSV file containing predictions.

    Returns:
    pd.DataFrame: A DataFrame containing the predictions for the specified UniProt ID.
    """
    uniprot_id = uniprot_id.upper()
    print(f"Fetching AlphaMissense predictions for {uniprot_id}")
    predictions = []

    with missense_tsv.open() as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == uniprot_id:
                prot_var = parts[1]
                predictions.append(
                    {
                        "protein_variant_from": prot_var[0],
                        "protein_variant_pos": int(prot_var[1:-1]),
                        "protein_variant_to": prot_var[-1],
                        "pathogenicity": float(parts[2]),
                        "classification": parts[3],
                    }
                )
    
    if not predictions:
        raise KeyError(f"No AlphaMissense predictions found for {uniprot_id}!")
    
    return pd.DataFrame(predictions)

if __name__ == "__main__":
    # Setup argparse to accept command-line arguments
    parser = argparse.ArgumentParser(description="Utility script for querying UniProt and downloading PDB files.")
    
    # Subparsers for different commands (UniProt query and PDB download)
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subparser for querying UniProt
    uniprot_parser = subparsers.add_parser("uniprot", help="Query UniProt for a gene's UniProt ID based on gene name and organism ID")
    uniprot_parser.add_argument('-n', '--gene-name', type=str, required=True, help="The gene name to search for (e.g., 'NAA10').")
    uniprot_parser.add_argument('-i', '--organism-id', type=int, required=True, help="The organism ID (e.g., 9606 for Homo sapiens).")

    # Subparser for downloading PDB files
    pdb_parser = subparsers.add_parser("pdb", help="Download a PDB file using its PDB ID")
    pdb_parser.add_argument('-p', '--pdb-id', type=str, required=True, help="The PDB ID of the protein structure (e.g., '6LID').")
    pdb_parser.add_argument('-o', '--output-dir', type=Path, required=True, help="Directory to save the downloaded PDB file.")

    args = parser.parse_args()

    # Handle the subcommands
    if args.command == "uniprot":
        uniprot_id = get_uniprot_id(args.gene_name, args.organism_id)
        if uniprot_id:
            print(f"UniProt ID for {args.gene_name} in organism {args.organism_id}: {uniprot_id}")
        else:
            print(f"No reviewed UniProt ID found for {args.gene_name} in organism {args.organism_id}.")
    
    elif args.command == "pdb":
        pdb_file_path = download_pdb_file(args.pdb_id, args.output_dir)
        if pdb_file_path:
            print(f"PDB file saved to: {pdb_file_path}")
        else:
            print("Failed to download the PDB file.")
