import gzip
import json
from pathlib import Path
from datetime import datetime
import numpy as np
import logging
import time
import Amissense.scripts.utils as utils

# Load configuration from config.json
config = utils.load_config()

def stream_tsv_file(tsv_file: Path, output_dir: Path):
    """
    Stream through the TSV file and generate JSON files for each uniprot_id.

    Parameters:
    tsv_file (Path): The path to the compressed TSV file.
    output_dir (Path): The path to the directory where JSON files will be saved.
    """
    utils.ensure_directory_exists(output_dir)

    # Initialize counters and timers
    file_count = 0
    start_time = time.time()

    with gzip.open(tsv_file, 'rt') as f:
        current_uniprot_id = None
        variants = {}
        am_pathogenicity_values = []
        am_class_counts = {}

        # Skip comment lines and look for the header
        for line in f:
            if not line.startswith('#'):
                header = line.strip().split('\t')
                break

        # Ensure that the header has the expected columns
        expected_columns = ['uniprot_id', 'protein_variant', 'am_pathogenicity', 'am_class']
        if header != expected_columns:
            raise ValueError(f"Unexpected header: {header}. Expected: {expected_columns}")

        # Read the file line by line
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')

            if len(parts) < 4:
                logging.warning(f"Skipping malformed line: {line.strip()}")
                continue

            uniprot_id = parts[0]
            protein_variant = parts[1]
            am_pathogenicity = float(parts[2])
            am_class = parts[3]

            if uniprot_id != current_uniprot_id:
                if current_uniprot_id is not None:
                    save_json_file(current_uniprot_id, variants, am_pathogenicity_values, am_class_counts, output_dir, tsv_file)
                    file_count += 1

                current_uniprot_id = uniprot_id
                variants = {protein_variant: {"am_pathogenicity": am_pathogenicity, "am_class": am_class}}
                am_pathogenicity_values = [am_pathogenicity]
                am_class_counts = {am_class: 1}
            else:
                variants[protein_variant] = {"am_pathogenicity": am_pathogenicity, "am_class": am_class}
                am_pathogenicity_values.append(am_pathogenicity)
                am_class_counts[am_class] = am_class_counts.get(am_class, 0) + 1

        if current_uniprot_id is not None:
            save_json_file(current_uniprot_id, variants, am_pathogenicity_values, am_class_counts, output_dir, tsv_file)
            file_count += 1

    logging.info(f"JSON file generation completed. Total files created: {file_count}. Time taken: {time.time() - start_time:.2f} seconds.")

def save_json_file(uniprot_id, variants, am_pathogenicity_values, am_class_counts, output_dir, tsv_file):
    """
    Save the list of variants and statistics to a JSON file for a specific uniprot_id.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    variants (dict): The dictionary of variants for this protein.
    am_pathogenicity_values (list): List of all am_pathogenicity values for calculating statistics.
    am_class_counts (dict): Dictionary counting the occurrences of each am_class.
    output_dir (Path): The directory where the JSON file will be saved.
    tsv_file (Path): The path to the source TSV file.
    """
    
    precision = config['defaults']['statistical_precision']
    am_pathogenicity_stats = {
        "average": round(np.mean(am_pathogenicity_values), precision),
        "min": round(np.min(am_pathogenicity_values), precision),
        "max": round(np.max(am_pathogenicity_values), precision),
        "median": round(np.median(am_pathogenicity_values), precision),
        "quantile_25": round(np.percentile(am_pathogenicity_values, 25), precision),
        "quantile_75": round(np.percentile(am_pathogenicity_values, 75), precision)
    }

    last_variant = list(variants.keys())[-1]
    protein_length = int(''.join(filter(str.isdigit, last_variant)))

    data = {
        "uniprot_id": uniprot_id,
        "variants": variants,
        "statistics": {
            "total_variants": len(variants),
            "am_class_counts": am_class_counts,
            "am_pathogenicity_stats": am_pathogenicity_stats,
            "protein_length": protein_length,
            "creation_date": datetime.now().isoformat(),
            "source_file": tsv_file.name
        }
    }

    json_file = output_dir / f"{uniprot_id}.AlphaMissense_aa_substitutions.json"
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)
    logging.info(f"Generated {json_file}")

if __name__ == "__main__":
    import argparse

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Generate JSON files from AlphaMissense TSV data.")
    parser.add_argument("tsv_file", type=Path, help="Path to the AlphaMissense_aa_substitutions.tsv.gz file.")
    parser.add_argument("output_dir", type=Path, help="Directory to store generated JSON files.")
    
    args = parser.parse_args()

    # Run the main function
    stream_tsv_file(args.tsv_file, args.output_dir)
