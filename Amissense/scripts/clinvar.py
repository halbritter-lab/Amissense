import requests
import pandas as pd
import re
from xml.etree import ElementTree
import time
import argparse
import Amissense.scripts.utils as utils

# Load configuration from config.json
config = utils.load_config()

# Dictionary to translate 3-letter amino acid codes to single-letter symbols
AA_DICT = config['amino_acid_codes']

def fetch_summary_ids(gene_id):
    url = f"{config['urls']['clinvar_esearch']}?db=clinvar&term={gene_id}[gene]&retmax=500"
    response = requests.get(url, timeout=5)
    response.raise_for_status()

    root = ElementTree.fromstring(response.content)
    id_list = root.find("IdList")
    return [id_elem.text for id_elem in id_list.findall("Id")]

def fetch_summaries(id_batch):
    ids = ",".join(id_batch)
    url = f"{config['urls']['clinvar_esummary']}?db=clinvar&id={ids}&retmode=json"
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()

def process_summaries(summaries):
    data = []
    result = summaries["result"]
    for uid in result["uids"]:
        variant_data = result[uid]

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
            refseq_number, gene_id, nucleotide_change = "", "", ""
            from_aa, location, to_aa = "", "", ""
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
    if batch_size is None:
        batch_size = config['defaults']['clinvar_batch_size']  # Use default from config.json
    print(f"Fetching ClinVar IDs for {gene_id}")
    try:
        ids = fetch_summary_ids(gene_id)
        print(f"Found {len(ids)} IDs.")
    except requests.HTTPError:
        print(f"Error fetching ClinVar IDs for {gene_id}")
        return None

    processed_summaries = []
    for index in range(0, len(ids), batch_size):
        try:
            time.sleep(1)  # Don't overwhelm ClinVar API
            summaries = fetch_summaries(ids[index : index + batch_size])
            processed_summaries.extend(process_summaries(summaries))
        except requests.HTTPError:
            print(f"Error fetching summaries for following IDs: {ids[index : index + batch_size]}")

    # Convert to dataframe and sort based on location
    output_df = pd.DataFrame(processed_summaries)
    output_df = output_df.sort_values("location")
    return output_df

def merge_missense_data(clinvar_data: pd.DataFrame, missense_data: pd.DataFrame) -> pd.DataFrame:
    # Merge and remove duplicate columns
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
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Fetch ClinVar data for a specific gene.")
    parser.add_argument("-g", "--gene-id", type=str, required=True, help="The Gene ID to fetch data for.")
    parser.add_argument("-b", "--batch-size", type=int, help="Number of IDs to fetch per batch.")
    
    args = parser.parse_args()

    # Fetch ClinVar data
    clinvar_data = fetch_clinvar_data(args.gene_id, batch_size=args.batch_size)

    # Display the fetched data if available
    if clinvar_data is not None:
        print(f"Fetched {clinvar_data.shape[0]} ClinVar entries for gene {args.gene_id}.")
        print(clinvar_data)
    else:
        print(f"No data fetched for gene {args.gene_id}.")
