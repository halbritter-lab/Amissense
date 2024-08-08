import requests
import pandas as pd
import re

# Dictionary to translate 3-letter amino acid codes to single-letter symbols
AA_DICT = {
    "Ala": "A",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y",
}


def fetch_summary_ids(gene_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_id}[gene]&retmax=500"
    response = requests.get(url)
    if response.status_code == 200:
        from xml.etree import ElementTree

        root = ElementTree.fromstring(response.content)
        id_list = root.find("IdList")
        return [id_elem.text for id_elem in id_list.findall("Id")]
    else:
        print(f"Error fetching IDs for {gene_id}: {response.status_code}")
        return []


def fetch_summaries(id_batch):
    ids = ",".join(id_batch)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={ids}&retmode=json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error fetching summaries: {response.status_code}")
        return None


def process_summaries(summaries):
    data = []
    result = summaries["result"]
    for uid in result["uids"]:
        variant_data = result[uid]

        if "missense variant" not in variant_data.get("molecular_consequence_list", []):
            continue

        germline_classification = variant_data.get("germline_classification", {}).get("description", "")

        variation_name = variant_data.get("variation_set", [{}])[0].get("variation_name", "")
        match = re.search(r"\(p\.([A-Za-z]+)(\d+)([A-Za-z]+)\)", variation_name)
        if match:
            from_aa, location, to_aa = match.groups()
            location = int(location)  # Convert location to integer
            from_aa_1letter = AA_DICT.get(from_aa, from_aa)  # Translate to 1-letter code
            to_aa_1letter = AA_DICT.get(to_aa, to_aa)  # Translate to 1-letter code
        else:
            from_aa, location, to_aa = "", "", ""
            from_aa_1letter, to_aa_1letter = "", ""

        data.append(
            {
                "germline_classification": germline_classification,
                "variation_name": variation_name,
                "location": location,
                "from": from_aa_1letter,
                "to": to_aa_1letter,
            }
        )

    return data


def fetch_clinvar_data(gene_id: str, batch_size=100) -> pd.DataFrame:
    ids = fetch_summary_ids(gene_id)
    if not ids:
        print("No IDs found.")
        return None

    print(f"Found {len(ids)} IDs.")

    processed_summaries = []
    for index in range(0, len(ids), batch_size):
        summaries = fetch_summaries(ids[index : index + batch_size])
        if summaries:
            processed_summaries.extend(process_summaries(summaries))

    # Convert to dataframe and sort based on location
    output_df = pd.DataFrame(processed_summaries)
    output_df = output_df.sort_values("location")

    return output_df


if __name__ == "__main__":
    GENE_ID = "SLC7A9"  # Replace with your desired protein ID
    df = fetch_clinvar_data(GENE_ID)
    print(df)
