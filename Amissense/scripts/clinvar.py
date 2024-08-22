import requests
import pandas as pd
import re
from xml.etree import ElementTree
import time
import logging
import Amissense.scripts.utils as utils

# Load configuration from config.json
try:
    config = utils.load_config()
except SystemExit:
    logging.critical("Unable to load configuration. Exiting.")
    raise

# Dictionary to translate 3-letter amino acid codes to single-letter symbols
AA_DICT = config['amino_acid_codes']

def fetch_summary_ids(gene_id):
    """
    Fetch ClinVar summary IDs for a given gene ID using E-Utilities API.
    
    Parameters:
    gene_id (str): The gene ID to search for.

    Returns:
    list: List of summary IDs.
    """
    url = f"{config['urls']['clinvar_esearch']}?db=clinvar&term={gene_id}[gene]&retmax=500"

    def make_request():
        response = requests.get(url, timeout=5)
        response.raise_for_status()
        return response

    try:
        response = utils.retry_request(make_request)
        root = ElementTree.fromstring(response.content)
        id_list = root.find("IdList")
        return [id_elem.text for id_elem in id_list.findall("Id")]
    except requests.RequestException as e:
        logging.error(f"Failed to fetch summary IDs for gene {gene_id}: {e}")
        raise

def fetch_summaries(id_batch):
    """
    Fetch ClinVar summaries for a batch of IDs using E-Utilities API.
    
    Parameters:
    id_batch (list): List of summary IDs to fetch.

    Returns:
    dict: JSON response from the E-Utilities API.
    """
    ids = ",".join(id_batch)
    url = f"{config['urls']['clinvar_esummary']}?db=clinvar&id={ids}&retmode=json"

    def make_request():
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response

    try:
        response = utils.retry_request(make_request)
        return response.json()
    except requests.RequestException as e:
        logging.error(f"Failed to fetch summaries for IDs {id_batch}: {e}")
        raise

def process_summaries(summaries):
    """
    Process the ClinVar summaries and extract relevant missense variant data.
    
    Parameters:
    summaries (dict): JSON response from the E-Utilities API.

    Returns:
    list: Processed ClinVar variant data.
    """
    data = []
    result = summaries.get("result", {})
    for uid in result.get("uids", []):
        variant_data = result.get(uid, {})

        if "missense variant" not in variant_data.get("molecular_consequence_list", []):
            continue

        germline_classification = variant_data.get("germline_classification", {}).get("description", "")

        variation_name = variant_data.get("variation_set", [{}])[0].get("variation_name", "")
        match = re.search(r"(NM_\d+\.\d+)\((\w+)\):c\.(\S+)\s+\(p\.([A-Za-z]+)(\d+)([A-Za-z]+)\)", variation_name)
        if match:
            refseq_number, gene_id, nucleotide_change, from_aa, location, to_aa = match.groups()
            location = int(location)  # Convert location to integer
            from_aa_1letter = AA_DICT.get(from_aa, from_aa)  # Translate to 1-letter code
            to_aa_1letter = AA_DICT.get(to_aa, to_aa)  # Translate to 1-letter code
        else:
            logging.warning(f"Skipping malformed variation name: {variation_name}")
            refseq_number, gene_id, nucleotide_change = "", "", ""
            from_aa_1letter, to_aa_1letter = "", ""

        data.append(
            {
                "germline_classification": germline_classification,
                "refseq_accession_number": refseq_number,
                "gene_id": gene_id,
                "nucleotide_change": nucleotide_change,
                "location": location,
                "from": from_aa_1letter,
                "to": to_aa_1letter,
            }
        )

    return data

def fetch_clinvar_data(gene_id: str, batch_size=None) -> pd.DataFrame:
    """
    Fetch ClinVar data for a specific gene ID in batches.

    Parameters:
    gene_id (str): The gene ID to fetch data for.
    batch_size (int, optional): The number of IDs to fetch per batch. Defaults to the value in config.

    Returns:
    pd.DataFrame: DataFrame containing processed ClinVar data.
    """
    if batch_size is None:
        batch_size = config['defaults']['clinvar_batch_size']  # Use default from config.json
    logging.info(f"Fetching ClinVar IDs for {gene_id}")
    try:
        ids = fetch_summary_ids(gene_id)
        logging.info(f"Found {len(ids)} IDs.")
    except requests.RequestException as e:
        logging.error(f"Error fetching ClinVar IDs for {gene_id}: {e}")
        return pd.DataFrame()  # Return an empty DataFrame on failure

    processed_summaries = []
    for index in range(0, len(ids), batch_size):
        try:
            time.sleep(1)  # Don't overwhelm ClinVar API
            summaries = fetch_summaries(ids[index: index + batch_size])
            processed_summaries.extend(process_summaries(summaries))
        except requests.RequestException as e:
            logging.error(f"Error fetching summaries for IDs {ids[index: index + batch_size]}: {e}")

    output_df = pd.DataFrame(processed_summaries)
    return output_df.sort_values("location") if not output_df.empty else output_df

def merge_missense_data(clinvar_data: pd.DataFrame, missense_data: pd.DataFrame) -> pd.DataFrame:
    """
    Merge ClinVar data with AlphaMissense data.

    Parameters:
    clinvar_data (pd.DataFrame): DataFrame containing ClinVar missense data.
    missense_data (pd.DataFrame): DataFrame containing AlphaMissense predictions.

    Returns:
    pd.DataFrame: Merged DataFrame.
    """
    merged_df = pd.merge(
        clinvar_data,
        missense_data,
        left_on=["from", "location", "to"],
        right_on=["protein_variant_from", "protein_variant_pos", "protein_variant_to"],
        how="left",
    )
    merged_df = merged_df.drop(["protein_variant_from", "protein_variant_pos", "protein_variant_to"], axis=1)

    # Rename columns
    merged_df = merged_df.rename(
        columns={
            "germline_classification": "ClinVar classification",
            "refseq_accession_number": "RefSeq accession number",
            "gene_id": "Gene",
            "nucleotide_change": "Nucleotide change (c.)",
            "location": "Amino acid change - location",
            "from": "Amino acid change - from",
            "to": "Amino acid change - to",
            "pathogenicity": "AM pathogenicity score",
            "classification": "AM classification",
        }
    )

    return merged_df

if __name__ == "__main__":
    import argparse

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Fetch ClinVar data for a specific gene.")
    parser.add_argument("-g", "--gene-id", type=str, required=True, help="The Gene ID to fetch data for.")
    parser.add_argument("-b", "--batch-size", type=int, help="Number of IDs to fetch per batch.")
    
    args = parser.parse_args()

    # Fetch ClinVar data
    try:
        clinvar_data = fetch_clinvar_data(args.gene_id, batch_size=args.batch_size)

        # Display the fetched data if available
        if not clinvar_data.empty:
            logging.info(f"Fetched {clinvar_data.shape[0]} ClinVar entries for gene {args.gene_id}.")
            logging.info(clinvar_data)
        else:
            logging.info(f"No data fetched for gene {args.gene_id}.")
    except Exception as e:
        logging.critical(f"Failed to fetch ClinVar data: {e}")
