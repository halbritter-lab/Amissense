

# Protein Pathogenicity Analysis Tool

This project provides tools to analyze and visualize the pathogenicity of protein variants based on AlphaMissense predictions.
It generates pathogenicity plots, modified PDB files, and integrates with ClinVar data to provide comprehensive insights into protein function and pathogenicity.


## Features

- Fetches and processes AlphaMissense predictions for specified proteins.
- Automatically fetches AlphaFold PDB files or utilizes user-provided experimental PDB files. #NEED TO AUTOMATE DOWNLOADING THE EXPERIMENTAL PDB FILES.
- Generates a modified PDB file with AlphaMissense pathogenicity scores.
    - Replaces the temperature factor (B-factor) with pathogenicity values, enabling visualization in molecular visualization tools.
- Creates a heatmap with the AlphaMissense predicted pathogenicity for all variants
- Produces a line graph visualizing:
    - Average predicted pathogenicity at each residue
    - AlphaFold confidence at each residue
    - Annotations for alpha helices and beta sheets
- Extracts all reported missense variants from ClinVar and integrates them into plots for comparisons with AlphaMissense predictions


## Requirements

- Python 3.7+
- Required Python packages:
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - requests
  - Biopython

You can install the required packages using:
    pip install pandas numpy matplotlib seaborn requests biopython


## Usage

1. Set the following constants in the script:
    - `UNIPROT_ID`: UniProt ID of the protein of interest
    - `GENE_ID`: Gene ID corresponding to the protein
    - `EXPERIMENTAL_PDB`: Path to an experimental PDB file (if available, otherwise leave as an empty string)

2. Run the script:
    python protein_pathogenicity_analysis.py #IS THIS TURE?


## Output

The script generates several output files in the `out/` directory:

- `{UNIPROT_ID}_pathogenicity_predictions.csv`: CSV file containing the AlphaMissense predictions
- `{UNIPROT_ID}_pathogenicity.pdb`: PDB file with pathogenicity scores within the B-factor
- `{UNIPROT_ID}_heatmap.png`: Heatmap of pathogenicity scores for all amino acid substitutions
- `{UNIPROT_ID}_line_graph.png`: Line graph showing mean pathogenicity and secondary structure information
- `{GENE_ID}_AM_clinvar.csv`: Combined data from AlphaMissense and ClinVar
- `{GENE_ID}_clinvar.png`: Graph comparing AlphaMissense predictions with ClinVar classifications


## Functions

- `fetch_predictions`: Downloads and extracts AlphaMissense predictions
- `fetch_predictions_for_id`: Retrieves predictions for a specific UniProt ID
- `fetch_alphafold_pdb`: Downloads AlphaFold PDB file for a given UniProt ID
- `extract_chain_id`: Extracts the relevant chain ID from a PDB file
- `generate_pathogenicity_pdb`: Creates a PDB file with pathogenicity scores
- `extract_secondary_structures`: Extracts secondary structure information from a PDB file
- `generate_graphs`: Creates heatmap and line graph visualizations
- `clinvar_graph`: Generates a graph comparing AlphaMissense predictions with ClinVar data

## Notes

- The script uses the `helpers.clinvar` module to fetch ClinVar data. Ensure this module is available in your project directory.
- DSSP is required for secondary structure prediction. Install it using:
    brew install brewsci/bio/dssp


## License

[Specify license here?]