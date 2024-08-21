import requests
import argparse

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

if __name__ == "__main__":
    # Setup argparse to accept command-line arguments
    parser = argparse.ArgumentParser(description="Query UniProt for a gene's UniProt ID based on gene name and organism ID.")
    parser.add_argument("gene_name", type=str, help="The gene name to search for (e.g., 'NAA10').")
    parser.add_argument("organism_id", type=int, help="The organism ID (e.g., 9606 for Homo sapiens).")
    
    args = parser.parse_args()

    # Call the function with the provided arguments
    uniprot_id = get_uniprot_id(args.gene_name, args.organism_id)

    if uniprot_id:
        print(f"UniProt ID for {args.gene_name} in organism {args.organism_id}: {uniprot_id}")
    else:
        print(f"No reviewed UniProt ID found for {args.gene_name} in organism {args.organism_id}.")
