import os
import requests
import gzip
import shutil
import json
import warnings
import pandas as pd
from pathlib import Path
import logging
from datetime import datetime
from xml.etree import ElementTree

def load_config(config_path: Path = None) -> dict:
    """
    Load the configuration from the specified JSON file or an environment variable.

    Parameters:
    config_path (Path): Path to the JSON configuration file.

    Returns:
    dict: The configuration settings as a dictionary.

    Raises:
    FileNotFoundError: If the config file is not found.
    json.JSONDecodeError: If the config file is malformed.
    """
    # Use environment variable if no config path is provided
    if config_path is None:
        config_path = Path(os.getenv("AMISSENSE_CONFIG_PATH", "Amissense/config.json"))

    try:
        with open(config_path, "r") as config_file:
            return json.load(config_file)
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {config_path}")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Malformed configuration file: {config_path} - {str(e)}")
        raise

# Load configuration at the module level with error handling
try:
    config = load_config()
except (FileNotFoundError, json.JSONDecodeError) as e:
    logging.critical(f"Failed to load configuration: {e}")
    raise SystemExit(1)

def ensure_directory_exists(directory: Path):
    """
    Ensure that the given directory exists. If not, create it.

    Parameters:
    directory (Path): The path to the directory.

    Raises:
    OSError: If the directory could not be created.
    """
    try:
        directory.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logging.error(f"Failed to create directory: {directory} - {str(e)}")
        raise

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

    Raises:
    requests.RequestException: If the network request fails.
    ValueError: If the response data is not in the expected format.
    """
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        if response.headers['Content-Type'] == 'application/json':
            return response.json()
        else:
            return ElementTree.fromstring(response.content)
    except requests.RequestException as e:
        logging.error(f"Failed to fetch data from {url}: {str(e)}")
        raise
    except ElementTree.ParseError as e:
        logging.error(f"Failed to parse XML from {url}: {str(e)}")
        raise ValueError(f"Invalid XML response from {url}")
    except json.JSONDecodeError as e:
        logging.error(f"Failed to parse JSON from {url}: {str(e)}")
        raise ValueError(f"Invalid JSON response from {url}")

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
    base_url = config["urls"]["uniprot_api"]
    query = f"(organism_id:{organism_id}) AND (gene:{gene_name})"
    fields = "accession,reviewed,id,gene_names,organism_name"
    params = {"query": query, "fields": fields, "format": "json"}

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        # Search for the "UniProtKB reviewed (Swiss-Prot)" entry
        for entry in data.get("results", []):
            if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
                return entry.get("primaryAccession")

        logging.warning(f"No reviewed (Swiss-Prot) entry found for gene: {gene_name}, organism: {organism_id}")
        return None

    except requests.RequestException as e:
        logging.error(f"Error querying UniProt API for gene: {gene_name}, organism: {organism_id} - {str(e)}")
        return None

def generate_output_directory(base_dir: Path, gene_id: str, uniprot_id: str) -> Path:
    """
    Generate a structured output directory name based on gene ID, UniProt ID, and the current date.

    Parameters:
    base_dir (Path): The base output directory.
    gene_id (str): The gene ID.
    uniprot_id (str): The UniProt ID.

    Returns:
    Path: The structured output directory path.
    """
    date_str = datetime.now().strftime("%Y-%m-%d")
    output_dir = base_dir / f"{gene_id}_{uniprot_id}_{date_str}"
    ensure_directory_exists(output_dir)
    return output_dir

def download_pdb(uniprot_id: str, output_dir: Path) -> Path:
    """
    Download a PDB file for a given UniProt ID. It first tries to fetch the file from AlphaFold; if unsuccessful, 
    it fetches it from the RCSB PDB website.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    output_dir (Path): The directory where the downloaded PDB file should be saved.

    Returns:
    Path: The path to the downloaded PDB file, or None if both downloads failed.
    """
    ensure_directory_exists(output_dir)
    date_str = datetime.now().strftime("%Y-%m-%d")
    alphafold_pdb_path = output_dir / f"{uniprot_id.upper()}_alphafold_{date_str}.pdb"

    try:
        if not alphafold_pdb_path.exists():
            alphafold_api = config['urls']['alphafold_api'].format(uniprot_id=uniprot_id.upper())
            logging.info(f"Downloading PDB file from AlphaFold API: {alphafold_api}")
            response = requests.get(alphafold_api, timeout=15)
            response.raise_for_status()
            data = response.json()
            pdb_url = data[0]["pdbUrl"]
            pdb_response = requests.get(pdb_url)
            pdb_response.raise_for_status()
            alphafold_pdb_path.write_bytes(pdb_response.content)
            logging.info(f"Download completed and saved to {alphafold_pdb_path}")
        return alphafold_pdb_path
    except requests.RequestException as e:
        logging.warning(f"Failed to download from AlphaFold: {str(e)}")

    # If AlphaFold download fails, try RCSB PDB
    pdb_id = uniprot_id[:4].upper()  # Assuming the first 4 characters represent the PDB ID
    pdb_file_path = output_dir / f"{pdb_id}_rcsb_{date_str}.pdb"
    pdb_url = config['urls']['pdb_download'].format(pdb_id=pdb_id)
    
    try:
        logging.info(f"Downloading PDB file from RCSB PDB: {pdb_url}")
        response = requests.get(pdb_url)
        response.raise_for_status()
        with open(pdb_file_path, "wb") as file:
            file.write(response.content)
        logging.info(f"Downloaded PDB file: {pdb_file_path}")
        return pdb_file_path
    except requests.RequestException as e:
        logging.error(f"Failed to download PDB file from RCSB PDB: {str(e)}")
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
    url = config["urls"]["alphamissense_predictions"]
    file_path = tmp_dir / Path(url).name
    tsv_path = file_path.with_suffix("")

    try:
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
    except requests.RequestException as e:
        logging.error(f"Failed to download AlphaMissense predictions: {str(e)}")
        raise
    except OSError as e:
        logging.error(f"Failed to extract AlphaMissense predictions: {str(e)}")
        raise

    return tsv_path

def get_predictions_from_json_for_uniprot_id(uniprot_id: str, json_dir: Path) -> pd.DataFrame:
    """
    Fetch AlphaMissense predictions for a specific UniProt ID from a JSON file.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    json_dir (Path): The path to the directory containing the AlphaMissense JSON files.

    Returns:
    pd.DataFrame: A DataFrame containing the predictions for the specified UniProt ID.

    Raises:
    FileNotFoundError: If the JSON file for the UniProt ID is not found.
    ValueError: If there is a mismatch between the UniProt ID in the file and the requested ID.
    """
    json_file = json_dir / f"{uniprot_id}.AlphaMissense_aa_substitutions.json"

    try:
        if not json_file.exists():
            raise FileNotFoundError(f"No JSON file found for {uniprot_id} at {json_file}")

        with open(json_file, "r") as file:
            data = json.load(file)

        if data["uniprot_id"].upper() != uniprot_id.upper():
            raise ValueError(f"UniProt ID mismatch: Expected {uniprot_id}, found {data['uniprot_id']}")

        predictions = [
            {
                "protein_variant_from": variant[0],
                "protein_variant_pos": int(variant[1:-1]),
                "protein_variant_to": variant[-1],
                "pathogenicity": variant_data["am_pathogenicity"],
                "classification": variant_data["am_class"],
            }
            for variant, variant_data in data["variants"].items()
        ]
        return pd.DataFrame(predictions)

    except (FileNotFoundError, ValueError) as e:
        logging.error(str(e))
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Failed to parse JSON file for {uniprot_id}: {str(e)}")
        raise

def get_predictions_from_static_api(uniprot_id: str) -> pd.DataFrame:
    """
    Fetch AlphaMissense predictions for a specific UniProt ID from the static JSON API.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.

    Returns:
    pd.DataFrame: A DataFrame containing the predictions for the specified UniProt ID.

    Raises:
    requests.RequestException: If the network request fails.
    """
    api_url = f"{config['urls']['static_json_api']}{uniprot_id}.AlphaMissense_aa_substitutions.json"
    
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        data = response.json()

        predictions = [
            {
                "protein_variant_from": variant[0],
                "protein_variant_pos": int(variant[1:-1]),
                "protein_variant_to": variant[-1],
                "pathogenicity": variant_data["am_pathogenicity"],
                "classification": variant_data["am_class"],
            }
            for variant, variant_data in data["variants"].items()
        ]
        return pd.DataFrame(predictions)
    except requests.RequestException as e:
        logging.error(f"Failed to fetch AlphaMissense predictions from the static API: {str(e)}")
        raise

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

    try:
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
            raise KeyError(f"No AlphaMissense predictions found for {uniprot_id}!")
        
        return pd.DataFrame(predictions)
    except (FileNotFoundError, KeyError) as e:
        logging.error(str(e))
        raise
