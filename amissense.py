import gzip
import shutil
import requests
import argparse
import pandas as pd
from pathlib import Path

import helpers.pdb
import helpers.graphs
import helpers.clinvar

# Constants -> REMOVE THIS LATER
# SLC7A9 __ P82251 __ 6lid.pdb
# SLC3A1 __ Q07837 __ 6lid.pdb
# CLDN16 __ Q9Y5I7
# python test.py P82251 SLC7A9 --output_dir out --experimental_pdb 6lid.pdb

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


def main(uniprot_id: str, gene_id: str, output_dir: Path, experimental_pdb: Path):
    ensure_directory_exist(output_dir)

    try:
        missense_tsv = download_predictions()
        predictions = fetch_predictions_for_id(uniprot_id, missense_tsv)
        predictions.to_csv(
            output_dir / f"{uniprot_id}_AM_pathogenicity_predictions.csv",
            index=False,
        )

        if experimental_pdb == Path():
            pdb_path = download_alphafold_pdb(uniprot_id)
        else:
            print("Using experimental PDB...")
            pdb_path = experimental_pdb
        chain_id = helpers.pdb.extract_chain_id(uniprot_id, pdb_path)
        helices, sheets = helpers.pdb.extract_secondary_structures(chain_id, pdb_path)
        pdb_confidences = (
            helpers.pdb.extract_positional_confidences(chain_id, pdb_path) if experimental_pdb == Path() else None
        )

        helpers.pdb.generate_pathogenicity_pdb(uniprot_id, predictions, pdb_path, output_dir)
        helpers.graphs.plot_predictions_heatmap(uniprot_id, predictions, output_dir)
        helpers.graphs.plot_predictions_line_graph(
            uniprot_id, predictions, output_dir, helices, sheets, pdb_confidences
        )

        clinvar_data = helpers.clinvar.fetch_clinvar_data(gene_id)
        clinvar_merged_data = helpers.clinvar.merge_missense_data(clinvar_data, predictions)
        clinvar_merged_data.to_csv(str(output_dir / f"{gene_id}_clinvar_AM.csv"), index=False)

        helpers.graphs.plot_clinvar_scatter(gene_id, predictions, clinvar_merged_data, output_dir)
        helpers.graphs.plot_clinvar_sankey(gene_id, predictions, clinvar_merged_data, output_dir)

    except requests.HTTPError as http_err:
        print(f"Unexpected error during download: {http_err}")
    except KeyError as key_err:
        print(key_err)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Process AlphaMissense predictions and generate visualizations.")
    parser.add_argument("uniprot_id", type=str, help="The UniProt ID of the protein.")
    parser.add_argument("gene_id", type=str, help="The Gene ID.")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="out",
        help="Directory to store output files.",
    )
    parser.add_argument(
        "--experimental_pdb",
        type=str,
        default="",
        help="Path to experimental PDB file.",
    )
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    main(
        args.uniprot_id,
        args.gene_id,
        Path(args.output_dir),
        Path(args.experimental_pdb),
    )
