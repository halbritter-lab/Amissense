import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio
import argparse
import logging
import Amissense.scripts.utils as utils
import matplotlib
from pathlib import Path
from typing import Optional
from plotly.subplots import make_subplots


matplotlib.use('Agg')  # Use a non-interactive backend

# Load configuration from config.json
config = utils.load_config()
plt.rc('font', size=20)

def plot_predictions_heatmap(uniprot_id: str, predictions: pd.DataFrame, out_dir: Path):
    logging.info("Generating heatmap graph...")

    # Structure predictions data
    unique_mutations = predictions["protein_variant_to"].unique()
    heatmap_data = np.zeros(
        shape=(
            len(unique_mutations),
            predictions["protein_variant_pos"].max() + 1,
        )
    )

    for row in predictions.itertuples(index=False, name="Row"):
        x = row.protein_variant_pos
        y = np.where(unique_mutations == row.protein_variant_to)[0][0]
        heatmap_data[y, x] = row.pathogenicity

    # Create graph
    fig1, ax1 = plt.subplots(figsize=tuple(config['defaults']['plot_figure_size']))

    colormap = config['defaults']['plot_colormap']
    heatmap = ax1.imshow(heatmap_data, aspect="auto", interpolation="none", cmap=colormap)
    yticks = np.arange(len(unique_mutations))

    ax1.set_yticks(yticks)
    ax1.set_yticklabels(unique_mutations)
    ax1.set_ylim(max(yticks)+0.5, min(yticks)-0.5)
    ax1.set_xlabel("Residue Sequence Number", fontsize=20)
    ax1.set_ylabel("Alternate Amino Acid", fontsize=20)
    ax1.set_title("Pathogenicity Heatmap", fontsize=24)

    # Add colorbar
    fig1.colorbar(heatmap, ax=ax1, shrink=0.6, label="AM Pathogenicity")


    # Save heatmap
    heatmap_path = out_dir / f"{uniprot_id}_heatmap.png"
    plt.savefig(heatmap_path, format="png", bbox_inches="tight")
    plt.close(fig1)

    logging.info(f"Heatmap stored as {heatmap_path}")


def plot_predictions_line_graph(
    uniprot_id: str,
    predictions: pd.DataFrame,
    out_dir: Path,
    alphafold_confidences: Optional[np.ndarray] = None,
):
    logging.info("Generating predictions line graph...")

    # Structure average positional pathogenicity
    predictions_grouped_means = predictions.groupby("protein_variant_pos")["pathogenicity"].mean()
    positional_means = predictions_grouped_means.reindex(
        range(0, predictions["protein_variant_pos"].max() + 1), fill_value=0
    ).to_numpy()

    # Generate and store line graph
    fig1, ax1 = plt.subplots(figsize=tuple(config['defaults']['plot_figure_size']))

    # Average pathogenicity
    ax1.plot(positional_means, label="Mean Pathogenicity", color="green")
    ax1.set_ylim(0, 1.01)
    ax1.margins(x=0, tight=True)
    ax1.set_xlabel("Residue Sequence Number", fontsize=20)
    ax1.set_ylabel("Mean Pathogenicity", fontsize=20)

    # Confidence
    if alphafold_confidences is not None:
        ax2 = ax1.twinx()
        ax2.plot(alphafold_confidences, label="AlphaFold Confidence", color="mediumvioletred")
        ax2.set_ylim(0, 1.01)
        ax2.set_ylabel("AlphaFold Confidence")
        ax2.margins(x=0, tight=True)

    # Secondary structures
    unique_structures = predictions.groupby('protein_variant_pos')['secondary_structure'].first().reset_index()
    unique_structures = unique_structures.sort_values('protein_variant_pos')
    helices = (unique_structures['secondary_structure'] == 'helix').astype(int).tolist()
    sheets = (unique_structures['secondary_structure'] == 'sheet').astype(int).tolist()
    secondary_structures_x = np.arange(len(helices))

    ax1.bar(
        secondary_structures_x,
        helices,
        label="Alpha helix",
        color="orange",
        alpha=0.25,
        width=1,
        align="center",
    )
    ax1.bar(
        secondary_structures_x,
        sheets,
        label="Beta sheet",
        color="purple",
        alpha=0.25,
        width=1,
        align="center",
    )

    # Labels and legend
    lines, labels = ax1.get_legend_handles_labels()
    if alphafold_confidences is not None:
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, bbox_to_anchor=(0.5, -0.09), ncol=4, fontsize=20)
    else:
        ax1.legend(lines, labels, bbox_to_anchor=(0.5, -0.09), ncol=4, fontsize=20)
    ax1.set_title("Mean Pathogenicity, AlphaFold Confidence, and Secondary Structures", fontsize=24)

    # Save line graph image
    line_graph_path = out_dir / f"{uniprot_id}_line_graph.png"
    plt.savefig(line_graph_path, format="png", bbox_inches="tight")
    plt.close(fig1)

    logging.info(f"Line graph stored as {line_graph_path}")


def plot_clinvar_scatter(
    gene_id: str, missense_data: pd.DataFrame, clinvar_missense_data: pd.DataFrame, out_dir: Path
):
    # Map classification to reduce number of labels
    clinvar_missense_data["clinvar_classification_mapped"] = clinvar_missense_data["ClinVar classification"].map(
        config['clinvar_classification_mapping']
    )

    # Structure average positional pathogenicity
    predictions_grouped_means = missense_data.groupby("protein_variant_pos")["pathogenicity"].mean()
    positional_means = predictions_grouped_means.reindex(
        range(0, missense_data["protein_variant_pos"].max() + 1), fill_value=0
    ).to_numpy()

    plt.figure(figsize=tuple(config['defaults']['scatter_figure_size']))

    # Line plot of positional pathogenicity means
    plt.plot(positional_means, label="Mean AlphaMissense Pathogenicity", color="green", alpha=0.5, linewidth=1)

    # Scatterplot
    sns.scatterplot(
        data=clinvar_missense_data,
        x="Amino acid change - location",
        y="AM pathogenicity score",
        hue="clinvar_classification_mapped",
        palette=config['defaults']['scatter_palette'],  # Use configurable color palette
        hue_order=config['defaults']['scatter_palette'].keys(),
        s=config['defaults']['scatter_marker_size'],  # Use configurable marker size
    )

    # Secondary structures
    unique_structures = missense_data.groupby('protein_variant_pos')['secondary_structure'].first().reset_index()
    unique_structures = unique_structures.sort_values('protein_variant_pos')
    helices = (unique_structures['secondary_structure'] == 'helix').astype(int).tolist()
    sheets = (unique_structures['secondary_structure'] == 'sheet').astype(int).tolist()
    secondary_structures_x = np.arange(len(helices))

    plt.bar(
        secondary_structures_x,
        helices,
        label="Alpha helix",
        color="orange",
        alpha=0.25,
        width=1,
        align="center",
    )

    plt.bar(
        secondary_structures_x,
        sheets,
        label="Beta sheet",
        color="purple",
        alpha=0.25,
        width=1,
        align="center",
    )

    # General plot settings
    plt.legend(bbox_to_anchor=(0, -0.15), loc="upper left", ncol=3, fontsize = 15)
    plt.title("AlphaMissense Predicted Pathogenicity with ClinVar Pathogenicity Classification", fontsize=24)
    plt.xlabel("Residue Sequence Number", fontsize=20)
    plt.ylabel("AlphaMissense Predicted Pathogenicity", fontsize=20)
    plt.ylim((0, 1.01))
    plt.xlim((0, len(positional_means)))


    clinvar_scatter_path = out_dir / f"{gene_id}_avgAM_clinvar.png"
    plt.savefig(clinvar_scatter_path, format="png", bbox_inches="tight")
    logging.info(f"ClinVar graph stored at {clinvar_scatter_path}")
    plt.close()


def plot_clinvar_sankey(gene_id: str, clinvar_missense_data: pd.DataFrame, out_dir: Path):
    # Map classification to reduce number of labels
    clinvar_missense_data["clinvar_classification_mapped"] = clinvar_missense_data["ClinVar classification"].map(
        config['clinvar_classification_mapping']
    )

    # Match AM labels to remapped ClinVar labels
    clinvar_missense_data["am_classification_mapped"] = clinvar_missense_data["AM classification"].map({
        "benign": "Likely benign",
        "ambiguous": "Uncertain significance",
        "pathogenic": "Likely pathogenic",
    })

    unival_df = (
        clinvar_missense_data.groupby(["clinvar_classification_mapped", "am_classification_mapped"])
        .size()
        .reset_index(name="count")
    )

    # Calculate percentages for ClinVar and AlphaMissense classifications
    total_count = unival_df["count"].sum()
    clinvar_percentages = unival_df.groupby("clinvar_classification_mapped")["count"].sum() / total_count * 100
    am_percentages = unival_df.groupby("am_classification_mapped")["count"].sum() / total_count * 100

    # Create custom labels with percentages for both sets of nodes
    clinvar_labels = [f"ClinVar: {label} ({clinvar_percentages[label]:.1f}%)" for label in clinvar_percentages.index]
    am_labels = [f"AM: {label} ({am_percentages[label]:.1f}%)" for label in am_percentages.index]
    node_labels = clinvar_labels + am_labels

    # Create mappings for custom names to indices
    input_mapping = {label: idx for idx, label in enumerate(clinvar_percentages.index)}
    output_mapping = {label: idx + len(input_mapping) for idx, label in enumerate(am_percentages.index)}

    # Update the DataFrame with new mappings
    unival_df["source"] = unival_df["clinvar_classification_mapped"].map(input_mapping)
    unival_df["target"] = unival_df["am_classification_mapped"].map(output_mapping)

    # Create Sankey diagram
    sankey_fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=config['defaults']['sankey_node_pad'],  # Configurable node padding
                    thickness=config['defaults']['sankey_node_thickness'],  # Configurable node thickness
                    line=dict(
                        color=config['defaults']['sankey_line_color'],  # Configurable line color
                        width=config['defaults']['sankey_line_width'],  # Configurable line width
                    ),
                    label=node_labels,
                    color=["lightblue"] * len(clinvar_labels) + ["lightgreen"] * len(am_labels),
                ),
                link=dict(
                    source=unival_df["source"],
                    target=unival_df["target"],
                    value=unival_df["count"],
                ),
            )
        ]
    )

    # Update layout
    sankey_fig.update_layout(
        title_text=f"{gene_id} ClinVar vs. AlphaMissense Classification",
        font_size=20,
        autosize=False,
        width=config['defaults']['sankey_plot_width'] * 1.2,  # Increase overall width
        height=config['defaults']['sankey_plot_height'],
        margin=dict(t=50, b=100, l=100, r=20)  # Adjust margins
    )

    # Save the combined figure as PNG
    sankey_path = out_dir / f"{gene_id}_sankey_diagram_with_stats.png"
    pio.write_image(sankey_fig, file=sankey_path, format="png", scale=2)  # Increase resolution
    logging.info(f"Sankey diagram with stats stored at {sankey_path}")

def main():
    parser = argparse.ArgumentParser(description="Generate visualizations for AlphaMissense predictions.")
    
    # Add arguments for the script
    parser.add_argument("-g", "--graph-type", type=str, required=True, choices=["heatmap", "line", "scatter", "sankey"], 
                        help="Type of graph to generate: heatmap, line, scatter, or sankey.")
    parser.add_argument("-u", "--uniprot-id", type=str, required=True, help="UniProt ID of the protein.")
    parser.add_argument("-o", "--out-dir", type=Path, required=True, help="Output directory for saving the graph.")
    parser.add_argument("-p", "--predictions", type=Path, required=True, help="Path to the predictions CSV file.")
    parser.add_argument("-c", "--clinvar-data", type=Path, help="Path to the ClinVar data CSV file (required for scatter and sankey).")
    parser.add_argument("--helices", type=Path, help="Path to the helices data file (required for line graph).")
    parser.add_argument("--sheets", type=Path, help="Path to the sheets data file (required for line graph).")
    parser.add_argument("--confidences", type=Path, help="Path to the AlphaFold confidences data file (optional for line graph).")
    
    args = parser.parse_args()

    # Load predictions data
    predictions = pd.read_csv(args.predictions)

    # Create output directory if it doesn't exist
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Call appropriate plotting function based on graph type
    if args.graph_type == "heatmap":
        plot_predictions_heatmap(args.uniprot_id, predictions, args.out_dir)
    elif args.graph_type == "line":
        if not args.helices or not args.sheets:
            logging.error("Error: Helices and sheets data are required for the line graph.")
            return
        helices = np.loadtxt(args.helices)
        sheets = np.loadtxt(args.sheets)
        confidences = np.loadtxt(args.confidences) if args.confidences else None
        plot_predictions_line_graph(args.uniprot_id, predictions, args.out_dir, helices, sheets, confidences)
    elif args.graph_type == "scatter" or args.graph_type == "sankey":
        if not args.clinvar_data:
            logging.error("Error: ClinVar data is required for scatter and sankey graphs.")
            return
        clinvar_data = pd.read_csv(args.clinvar_data)
        if args.graph_type == "scatter":
            plot_clinvar_scatter(args.uniprot_id, predictions, clinvar_data, args.out_dir)
        else:
            plot_clinvar_sankey(args.uniprot_id, clinvar_data, args.out_dir)
    else:
        logging.error(f"Unknown graph type: {args.graph_type}")


if __name__ == "__main__":
    main()