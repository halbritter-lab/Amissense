import argparse
import gzip
import shutil
import requests
import pandas as pd
from pathlib import Path

import Amissense.scripts.pdb as pdb_module
import Amissense.scripts.graphs as graphs_module
import Amissense.scripts.clinvar as clinvar_module

TMP_DIR = Path("tmp")


def ensure_directory_exist(dir: Path):
    """Create directories if they don't exist."""
    dir.mkdir(parents=True, exist_ok=True)


def download_predictions() -> Path:
    """Download and extract AlphaMissense predictions."""
    ensure_directory_exist(TMP_DIR)

    url = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
    file_path = TMP_DIR / Path(url).name
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


def fetch_predictions_for_id(uniprot_id: str, missense_tsv: Path) -> pd.DataFrame:
    """Fetch AlphaMissense predictions for a specific UniProt ID."""
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


def download_alphafold_pdb(uniprot_id: str) -> Path:
    """Fetch AlphaFold PDB file for a given UniProt ID."""
    download_dir = TMP_DIR / "pdb_files"
    ensure_directory_exist(download_dir)
    alphafold_pdb_path = download_dir / f"{uniprot_id.upper()}_alphafold.pdb"

    if not alphafold_pdb_path.exists():
        alphafold_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"

        print(f"Downloading PDB file from {alphafold_api}")
        response = requests.get(alphafold_api, timeout=15)
        response.raise_for_status()
        data = response.json()

        pdb_url = data[0]["pdbUrl"]
        pdb_response = requests.get(pdb_url)
        pdb_response.raise_for_status()

        alphafold_pdb_path.write_bytes(pdb_response.content)
        print(f"Download completed and saved to {alphafold_pdb_path}")

    return alphafold_pdb_path


def run_pipeline(uniprot_id: str, gene_id: str, output_dir: Path, experimental_pdb: Path):
    """Main function to run the pipeline."""
    ensure_directory_exist(output_dir)

    try:
        # Download and fetch AlphaMissense predictions
        missense_tsv = download_predictions()
        predictions = fetch_predictions_for_id(uniprot_id, missense_tsv)
        predictions.to_csv(
            output_dir / f"{uniprot_id}_AM_pathogenicity_predictions.csv",
            index=False,
        )

        # Use experimental PDB if provided, otherwise download from AlphaFold
        if experimental_pdb == Path():
            pdb_path = download_alphafold_pdb(uniprot_id)
        else:
            print("Using experimental PDB...")
            pdb_path = experimental_pdb

        # Extract PDB details
        chain_id = pdb_module.extract_chain_id(uniprot_id, pdb_path)
        helices, sheets = pdb_module.extract_secondary_structures(chain_id, pdb_path)
        pdb_confidences = (
            pdb_module.extract_positional_confidences(chain_id, pdb_path) if experimental_pdb == Path() else None
        )

        # Generate PDB with pathogenicity and create visualizations
        pdb_module.generate_pathogenicity_pdb(uniprot_id, predictions, pdb_path, output_dir)
        graphs_module.plot_predictions_heatmap(uniprot_id, predictions, output_dir)
        graphs_module.plot_predictions_line_graph(
            uniprot_id, predictions, output_dir, helices, sheets, pdb_confidences
        )

        # Fetch and merge ClinVar data, then create visualizations
        clinvar_data = clinvar_module.fetch_clinvar_data(gene_id)
        clinvar_merged_data = clinvar_module.merge_missense_data(clinvar_data, predictions)
        clinvar_merged_data.to_csv(str(output_dir / f"{gene_id}_clinvar_AM.csv"), index=False)

        graphs_module.plot_clinvar_scatter(gene_id, predictions, clinvar_merged_data, output_dir)
        graphs_module.plot_clinvar_sankey(gene_id, predictions, clinvar_merged_data, output_dir)

    except requests.HTTPError as http_err:
        print(f"Unexpected error during download: {http_err}")
    except KeyError as key_err:
        print(key_err)


if __name__ == "__main__":
    # Setup argument parsing
    parser = argparse.ArgumentParser(description="Run the AlphaMissense data processing pipeline.")
    parser.add_argument('-u', '--uniprot-id', type=str, required=True, help="The UniProt ID of the protein.")
    parser.add_argument('-g', '--gene-id', type=str, required=True, help="The Gene ID.")
    parser.add_argument('-o', '--output-dir', type=str, default="out", help="Directory to store output files.")
    parser.add_argument('-e', '--experimental-pdb', type=str, default="", help="Path to experimental PDB file.")
    
    args = parser.parse_args()

    # Call the main pipeline function
    run_pipeline(args.uniprot_id, args.gene_id, Path(args.output_dir), Path(args.experimental_pdb))
