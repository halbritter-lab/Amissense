# Protein Pathogenicity Analysis Tool

This project offers an accessible script to analyze and visualize AlphaMissense-predicted pathogenicity scores. 
It integrates with ClinVar and AlphaFold to promote a more comprehensive interpretation of the pathogenicity predictions.

## Features

- Fetches and processes AlphaMissense predictions for a specified protein (based on UniProt ID)
- Automatically fetches AlphaFold PDB files or utilizes user-provided experimental PDB files to generate a modified PDB file with AlphaMissense pathogenicity scores
  - Replaces the temperature factor (B-factor) with pathogenicity values, enabling visualization in molecular visualization tools
- Creates a heatmap with the AlphaMissense predicted pathogenicity scores for all amino acid residue substitutions
- Extracts reporter missense variants from ClinVar
  - Compares the pathogenicity classification of all missense ClinVar variants to the AlphaMissense classifications (exported as `{GENE_ID}_clinvar_AM.csv`)
- Produces a line graph visualizing:
  - Average AlphaMissense predicted pathogenicity at each amino acid residue
  - AlphaFold per-residue model confidence score (pLDDT)
  - Secondary structure annotations for alpha helices and beta sheets
- Produces another line graph visualizing:
  - Average AlphaMissense predicted pathogenicity at each amino acid residue
  - Missense variants from ClinVar along with their associated classification (simplified for plot simplicity)

## Requirements

- Python 3.11 or higher

You can install the required packages using: 
```bash
pip install -r requirements.txt
```

Alternatively, you can create a conda environment with the required packages using:
```bash
conda create --name amissense python requests pandas seaborn matplotlib plotly numpy biopython conda-forge::python-kaleido salilab::dssp
conda activate amissense
```

## Installation

With the `setup.py` file, you can now install the package using the following:

```bash
pip install .
```

This will install the `amissense` package and make the `amissense` command-line interface (CLI) available.

## Usage

To run the main pipeline, use the following command:

```bash
amissense pipeline -u UNIPROT_ID -g GENE_ID [-o OUTPUT_DIR] [-e EXPERIMENTAL_PDB]
```

Arguments:
- `UNIPROT_ID`: The UniProt ID of the protein you want to analyze. This is a required positional argument.
- `GENE_ID`: The gene ID associated with the protein. This is a required positional argument.
- `OUTPUT_DIR`: The directory to store output files (default is `out/`).
- `EXPERIMENTAL_PDB`: The experimental PDB file for the protein (optional). If not provided, the script will use the AlphaFold predicted structure.

### Example Command

```bash
amissense pipeline -u P12345 -g BRCA1 -o output_dir -e experimental.pdb
```

### Utility Commands

You can also use utility commands for additional operations:

- **Download AlphaMissense predictions:**
```bash
amissense utils download-predictions -t /path/to/tmp_dir
```

- **Download a PDB file:**
```bash
amissense utils download-pdb -p 6LID -o /path/to/output_dir
```

- **Query UniProt for a gene's UniProt ID:**
```bash
amissense utils uniprot-query -n GENE_NAME -i ORGANISM_ID
```

### JSON Generation

To generate JSON files from AlphaMissense TSV data:
```bash
amissense generate-json /path/to/AlphaMissense_aa_substitutions.tsv.gz /path/to/output_dir
```

## Output

The script generates several output files in the specified output directory:

- `{UNIPROT_ID}_heatmap.png`: Heatmap of the AlphaMissense pathogenicity scores for each amino acid substitution
- `{UNIPROT_ID}_line_graph.png`: Line graph showing the AlphaMissense mean pathogenicity, AlphaFold per-residue model confidence score (pLDDT), and secondary structure annotations
- `{UNIPROT_ID}_pathogenicity.pdb`: PDB file with pathogenicity scores within the B-factor
- `{UNIPROT_ID}_AM_pathogenicity_predictions.csv`: CSV file containing the AlphaMissense predictions
- `{GENE_ID}_clinvar_AM.csv`: CSV file with ClinVar and AlphaMissense data combined
- `{GENE_ID}_avgAM_clinvar.png`: Line graph showing the AlphaMissense mean pathogenicity and extracted ClinVar missense variants with their classification
- `{GENE_ID}_sankey_diagram.png`: Sankey diagram depicting the flow quantity between the AlphaMissense variant pathogenicity classification and the ClinVar variant classifications
