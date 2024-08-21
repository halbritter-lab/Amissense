import gzip
import json
from pathlib import Path
from datetime import datetime
import numpy as np
import logging
import time

def ensure_directory_exists(directory: Path):
    """
    Ensure that the given directory exists. If not, create it.

    Parameters:
    directory (Path): The path to the directory.
    """
    directory.mkdir(parents=True, exist_ok=True)

def stream_tsv_file(tsv_file: Path, output_dir: Path):
    """
    Stream through the TSV file and generate JSON files for each uniprot_id.

    Parameters:
    tsv_file (Path): The path to the compressed TSV file.
    output_dir (Path): The path to the directory where JSON files will be saved.
    """
    ensure_directory_exists(output_dir)

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
            # Skip comment lines or empty lines
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')

            # Ensure the line has enough parts (columns)
            if len(parts) < 4:
                logging.warning(f"Skipping malformed line: {line.strip()}")
                continue

            uniprot_id = parts[0]
            protein_variant = parts[1]
            am_pathogenicity = float(parts[2])
            am_class = parts[3]

            # If we encounter a new uniprot_id, write the current data to a JSON file
            if uniprot_id != current_uniprot_id:
                if current_uniprot_id is not None:
                    file_start_time = time.time()  # Start timing for this file
                    save_json_file(current_uniprot_id, variants, am_pathogenicity_values, am_class_counts, output_dir, tsv_file)
                    file_count += 1
                    logging.info(f"Generated {file_count} files so far.")
                    logging.info(f"Time taken for last file: {time.time() - file_start_time:.2f} seconds")

                # Start collecting data for the new uniprot_id
                current_uniprot_id = uniprot_id
                variants = {protein_variant: {"am_pathogenicity": am_pathogenicity, "am_class": am_class}}
                am_pathogenicity_values = [am_pathogenicity]
                am_class_counts = {am_class: 1}
            else:
                # Add the variant to the current uniprot_id's list
                variants[protein_variant] = {"am_pathogenicity": am_pathogenicity, "am_class": am_class}
                am_pathogenicity_values.append(am_pathogenicity)

                # Count occurrences of am_class
                if am_class in am_class_counts:
                    am_class_counts[am_class] += 1
                else:
                    am_class_counts[am_class] = 1

        # Write the last uniprot_id's data to a JSON file
        if current_uniprot_id is not None:
            file_start_time = time.time()  # Start timing for this file
            save_json_file(current_uniprot_id, variants, am_pathogenicity_values, am_class_counts, output_dir, tsv_file)
            file_count += 1
            logging.info(f"Generated {file_count} files so far.")
            logging.info(f"Time taken for last file: {time.time() - file_start_time:.2f} seconds")

    # Log the total time taken for all files
    total_time = time.time() - start_time
    logging.info(f"JSON file generation completed. Total files created: {file_count}. Time taken: {total_time:.2f} seconds.")

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
    
    # Calculate statistics with rounding
    total_variants = len(variants)
    am_pathogenicity_stats = {
        "average": round(np.mean(am_pathogenicity_values), 4),
        "min": round(np.min(am_pathogenicity_values), 4),
        "max": round(np.max(am_pathogenicity_values), 4),
        "median": round(np.median(am_pathogenicity_values), 4),
        "quantile_25": round(np.percentile(am_pathogenicity_values, 25), 4),
        "quantile_75": round(np.percentile(am_pathogenicity_values, 75), 4)
    }

    # Extract protein length from the last variant's position (e.g., Q1630S -> 1630)
    last_variant = list(variants.keys())[-1]
    protein_length = int(''.join(filter(str.isdigit, last_variant)))

    # Date of creation
    creation_date = datetime.now().isoformat()

    # Extract source file name
    source_file = tsv_file.name

    data = {
        "uniprot_id": uniprot_id,
        "variants": variants,
        "statistics": {
            "total_variants": total_variants,
            "am_class_counts": am_class_counts,
            "am_pathogenicity_stats": am_pathogenicity_stats,
            "protein_length": protein_length,
            "creation_date": creation_date,
            "source_file": source_file
        }
    }

    # Define the JSON file path
    json_file = output_dir / f"{uniprot_id}.AlphaMissense_aa_substitutions.json"

    # Write the JSON file
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)
    logging.info(f"Generated {json_file}")

def main(tsv_file: Path, output_dir: Path):
    """
    Main function to parse the TSV file and generate JSON files in a streaming manner.

    Parameters:
    tsv_file (Path): The path to the compressed TSV file.
    output_dir (Path): The path to the directory where JSON files will be saved.
    """
    logging.info(f"Processing {tsv_file}...")
    stream_tsv_file(tsv_file, output_dir)

if __name__ == "__main__":
    import argparse

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Generate JSON files from AlphaMissense TSV data.")
    parser.add_argument("tsv_file", type=Path, help="Path to the AlphaMissense_aa_substitutions.tsv.gz file.")
    parser.add_argument("output_dir", type=Path, help="Directory to store generated JSON files.")
    
    args = parser.parse_args()

    # Run the main function
    main(args.tsv_file, args.output_dir)
