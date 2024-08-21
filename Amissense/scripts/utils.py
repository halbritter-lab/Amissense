import requests
import gzip
import shutil
import json
import warnings
import pandas as pd
from pathlib import Path
import logging
from xml.etree import ElementTree

def load_config(config_path: Path = Path("Amissense/config.json")) -> dict:
    """
    Load the configuration from the specified JSON file.

    Parameters:
    config_path (Path): Path to the JSON configuration file.

    Returns:
    dict: The configuration settings as a dictionary.
    """
    with open(config_path, "r") as config_file:
        return json.load(config_file)

# Load configuration at the module level
config = load_config()

def ensure_directory_exists(directory: Path):
    """
    Ensure that the given directory exists. If not, create it.

    Parameters:
    directory (Path): The path to the directory.
    """
    directory.mkdir(parents=True, exist_ok=True)

def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Set up logging configuration.

    Parameters:
    log_level: Logging level (default: INFO)
    log_file: Optional file path to log to a file
    """
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        filename=log_file,
        filemode='a'
    )
    if not log_file:
        logging.getLogger().addHandler(logging.StreamHandler())

def fetch_data_from_url(url: str, timeout: int = 5) -> dict:
    """
    Fetch data from a given URL and return it as a parsed dictionary.

    Parameters:
    url (str): The URL to fetch data from.
    timeout (int): Timeout for the request in seconds.

    Returns:
    dict: Parsed JSON or XML data from the response.
    """
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.json() if response.headers['Content-Type'] == 'application/json' else ElementTree.fromstring(response.content)

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
    base_url = config["apis"]["uniprot"]["url"]
    query = f"(organism_id:{organism_id}) AND (gene:{gene_name})"
    fields = config["apis"]["uniprot"]["fields"]
    params = {"query": query, "fields": fields, "format": "json"}

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        for entry in data.get("results", []):
            if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
                return entry.get("primaryAccession")
    except requests.RequestException as e:
        logging.error(f"Error querying UniProt API: {e}")

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
    ensure_directory_exists(output_dir)
    pdb_url = f"{config['apis']['pdb']['url']}/{pdb_id}.pdb"
    
    try:
        response = requests.get(pdb_url)
        response.raise_for_status()
        pdb_file_path = output_dir / f"{pdb_id}.pdb"
        with open(pdb_file_path, "wb") as file:
            file.write(response.content)
        logging.info(f"Downloaded PDB file: {pdb_file_path}")
        return pdb_file_path
    except requests.RequestException as e:
        logging.error(f"Error downloading PDB file {pdb_id}: {e}")
        return None

def download_and_extract_alphamissense_predictions(tmp_dir: Path) -> Path:
    """
    Download and extract AlphaMissense predictions.

    Parameters:
    tmp_dir (Path): The temporary directory where the prediction file will be downloaded and extracted.

    Returns:
    Path: The path to the extracted TSV file.
    """
    ensure_directory_exists(tmp_dir)
    url = config["apis"]["alphamissense"]["url"]
    file_path = tmp_dir / Path(url).name
    tsv_path = file_path.with_suffix("")

    if not tsv_path.exists():
        if not file_path.exists():
            logging.info(f"Downloading AlphaMissense predictions from {url}")
            response = requests.get(url)
            response.raise_for_status()
            file_path.write_bytes(response.content)
            logging.info("Download completed!")
        
        logging.info("Extracting predictions")
        with gzip.open(file_path, "rb") as f_in, tsv_path.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logging.info("Extraction completed!")
    
    return tsv_path

def get_predictions_from_json_for_uniprot_id(uniprot_id: str, json_dir: Path) -> pd.DataFrame:
    """
    Fetch AlphaMissense predictions for a specific UniProt ID from a JSON file.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    json_dir (Path): The path to the directory containing the AlphaMissense JSON files.

    Returns:
    pd.DataFrame: A DataFrame containing the predictions for the specified UniProt ID.
    """
    json_file = json_dir / f"{uniprot_id}.AlphaMissense_aa_substitutions.json"

    if not json_file.exists():
        logging.error(f"No JSON file found for {uniprot_id} at {json_file}")
        raise FileNotFoundError(f"No JSON file found for {uniprot_id} at {json_file}")

    with open(json_file, "r") as file:
        data = json.load(file)

    if data["uniprot_id"].upper() != uniprot_id.upper():
        logging.error(f"UniProt ID mismatch: Expected {uniprot_id}, found {data['uniprot_id']}")
        raise ValueError(f"UniProt ID mismatch: Expected {uniprot_id}, found {data['uniprot_id']}")

    predictions = []
    for variant, variant_data in data["variants"].items():
        predictions.append(
            {
                "protein_variant_from": variant[0],
                "protein_variant_pos": int(variant[1:-1]),
                "protein_variant_to": variant[-1],
                "pathogenicity": variant_data["am_pathogenicity"],
                "classification": variant_data["am_class"],
            }
        )
    return pd.DataFrame(predictions)

def get_predictions_from_am_tsv_for_uniprot_id(uniprot_id: str, missense_tsv: Path) -> pd.DataFrame:
    """
    Deprecated: Fetch AlphaMissense predictions for a specific UniProt ID from a TSV file.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    missense_tsv (Path): The path to the AlphaMissense TSV file containing predictions.

    Returns:
    pd.DataFrame: A DataFrame containing the predictions for the specified UniProt ID.
    """
    warnings.warn(
        "get_predictions_from_am_tsv_for_uniprot_id is deprecated. Use get_predictions_from_json_for_uniprot_id instead.",
        DeprecationWarning
    )

    uniprot_id = uniprot_id.upper()
    logging.info(f"Fetching AlphaMissense predictions for {uniprot_id}")
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
        logging.error(f"No AlphaMissense predictions found for {uniprot_id}!")
        raise KeyError(f"No AlphaMissense predictions found for {uniprot_id}!")
    
    return pd.DataFrame(predictions)
