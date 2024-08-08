import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from typing import Optional
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio


def plot_predictions_heatmap(uniprot_id: str, predictions: pd.DataFrame, out_dir: Path):
    print("Generating heatmap graph...")

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
    fig1, ax1 = plt.subplots(figsize=(30, 10))

    colormap = "coolwarm"
    heatmap = ax1.imshow(heatmap_data, aspect="auto", interpolation="none", cmap=colormap)
    yticks = np.arange(len(unique_mutations)) + 0.5

    ax1.set_yticks(yticks)
    ax1.set_yticklabels(unique_mutations)
    ax1.set_ylim(20, 0)  # Adjust as needed
    ax1.set_xlabel("Residue Sequence Number")
    ax1.set_ylabel("Alternate Amino Acid")
    ax1.set_title("Pathogenicity Heatmap")

    # Add colorbar
    fig1.colorbar(heatmap, ax=ax1, shrink=0.6, label="AM Pathogenicity")

    # Save heatmap
    heatmap_path = out_dir / f"{uniprot_id}_heatmap.png"
    plt.savefig(heatmap_path, format="png", bbox_inches="tight")
    plt.close(fig1)

    print(f"Heatmap stored as {heatmap_path}")


def plot_predictions_line_graph(
    uniprot_id: str,
    predictions: pd.DataFrame,
    out_dir: Path,
    helices: np.ndarray,
    sheets: np.ndarray,
    alphafold_confidences: Optional[np.ndarray] = None,
):
    print("Generating predictions line graph...")

    # Structure average positional pathogenicity
    predictions_grouped_means = predictions.groupby("protein_variant_pos")["pathogenicity"].mean()
    positional_means = predictions_grouped_means.reindex(
        range(0, predictions["protein_variant_pos"].max() + 1), fill_value=0
    ).to_numpy()

    # Generate and store line graph
    fig1, ax1 = plt.subplots(figsize=(30, 10))

    # Average pathogenicity
    ax1.plot(positional_means, label="Mean Pathogenicity", color="green")
    ax1.set_ylim(0, 1.1)
    ax1.set_xlabel("Residue Sequence Number")
    ax1.set_ylabel("Mean Pathogenicity")
    ax1.margins(x=0, tight=True)

    # Confidence
    if alphafold_confidences is not None:
        ax2 = ax1.twinx()
        ax2.plot(alphafold_confidences, label="AlphaFold Confidence", color="mediumvioletred")
        ax2.set_ylim(0, 1.1)
        ax2.set_ylabel("AlphaFold Confidence")
        ax2.margins(x=0, tight=True)

    # Secondary structures
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
        ax1.legend(lines + lines2, labels + labels2, loc="upper left")
    else:
        ax1.legend(lines, labels, loc="upper left")
    ax1.set_title("Mean Pathogenicity, AlphaFold Confidence, and Secondary Structures")

    # Save line graph image
    line_graph_path = out_dir / f"{uniprot_id}_line_graph.png"
    plt.savefig(line_graph_path, format="png", bbox_inches="tight")
    plt.close(fig1)

    print(f"Line graph stored as {line_graph_path}")


def _germline_classification_mapping() -> dict:
    return {
        "Pathogenic": "Likely pathogenic",
        "Pathogenic/Likely pathogenic": "Likely pathogenic",
        "Likely pathogenic": "Likely pathogenic",
        "Benign": "Likely benign",
        "Benign/Likely benign": "Likely benign",
        "Likely benign": "Likely benign",
        "Uncertain significance": "Uncertain significance",
        "Conflicting classifications of pathogenicity": "Conflicting significance",
    }


def _pathogenicity_mapping(value) -> str:
    if 0 <= value <= 0.34:
        return "Likely benign"
    elif 0.35 <= value <= 0.564:
        return "Uncertain significance"
    elif 0.565 <= value <= 1:
        return "Likely pathogenic"
    else:
        return "Unknown"


def plot_clinvar_scatter(gene_id: str, predictions: pd.DataFrame, clinvar_data: pd.DataFrame, out_dir: Path):
    # Merge AM_pathogenicity scores to clinvar df
    merged_data = pd.merge(
        clinvar_data,
        predictions,
        left_on=["from", "location", "to"],
        right_on=["protein_variant_from", "protein_variant_pos", "protein_variant_to"],
        how="left",
    )

    # Apply classification mapping
    merged_data["classification_mapped"] = merged_data["germline_classification"].map(
        _germline_classification_mapping()
    )

    # Save the merged data
    merged_data.to_csv(str(out_dir / f"{gene_id}_AM_clinvar_merged.csv"), index=False)

    # Structure average positional pathogenicity
    predictions_grouped_means = predictions.groupby("protein_variant_pos")["pathogenicity"].mean()
    positional_means = predictions_grouped_means.reindex(
        range(0, predictions["protein_variant_pos"].max() + 1), fill_value=0
    ).to_numpy()

    plt.figure(figsize=(20, 6))

    # Scatterplot with custom color palette
    pathogenicity_palette = {
        "Likely pathogenic": "orangered",
        "Likely benign": "royalblue",
        "Uncertain significance": "orange",
        "Conflicting significance": "gray",
    }

    sns.scatterplot(
        data=merged_data,
        x="location",
        y="pathogenicity",
        hue="classification_mapped",
        palette=pathogenicity_palette,
        hue_order=pathogenicity_palette.keys(),
        s=100,
    )

    # Line plot of positional pathogenicity means
    plt.plot(positional_means, label="Mean AlphaMissense Pathogenicity", color="green", alpha=0.5, linewidth=1)

    # General plot settings
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    plt.title("AlphaMissense Predicted Pathogenicity with Clinvar Pathogenicity Classification")
    plt.xlabel("Residue Sequence Number")
    plt.ylabel("AlphaMissense Predicted Pathogenicity")

    clinvar_scatter_path = out_dir / f"{gene_id}_avgAM_clinvar.png"
    plt.savefig(clinvar_scatter_path, format="png", bbox_inches="tight")
    print(f"Clinvar graph stored at {clinvar_scatter_path}")
    plt.close()


def plot_clinvar_sankey(gene_id: str, predictions: pd.DataFrame, clinvar_data: pd.DataFrame, out_dir: Path):
    # Merge AM_pathogenicity scores to clinvar df
    merged_data = pd.merge(
        clinvar_data,
        predictions,
        left_on=["from", "location", "to"],
        right_on=["protein_variant_from", "protein_variant_pos", "protein_variant_to"],
        how="left",
    )

    # Create the clinvar_amissense DataFrame
    clinvar_amissense = merged_data[["germline_classification", "pathogenicity"]].copy()

    # Apply the classification mappings
    clinvar_amissense["classification_mapped"] = clinvar_amissense["germline_classification"].map(
        _germline_classification_mapping()
    )
    clinvar_amissense["pathogenicity_mapped"] = clinvar_amissense["pathogenicity"].apply(_pathogenicity_mapping)

    # Append new column to check if classification_mapped equals pathogenicity_mapped
    clinvar_amissense["match"] = clinvar_amissense.apply(
        lambda row: "yes" if row["classification_mapped"] == row["pathogenicity_mapped"] else "no", axis=1
    )

    # Save the clinvar_amissense DataFrame as a CSV file
    clinvar_amissense.to_csv(str(out_dir / f"{gene_id}_AM_clinvar_matched.csv"), index=False)

    # Calculate and print the percentages of 'yes' and 'no'
    match_counts = clinvar_amissense["match"].value_counts(normalize=True) * 100
    print("Some info on what's printed below...")  # TODO
    print(f"Percentage of 'yes': {match_counts.get('yes', 0):.2f}%")
    print(f"Percentage of 'no': {match_counts.get('no', 0):.2f}%")

    unival_df = (
        clinvar_amissense.groupby(["classification_mapped", "pathogenicity_mapped"]).size().reset_index(name="count")
    )

    # Calculate percentages for ClinVar and AlphaMissense classifications
    total_count = unival_df["count"].sum()
    clinvar_percentages = unival_df.groupby("classification_mapped")["count"].sum() / total_count * 100
    am_percentages = unival_df.groupby("pathogenicity_mapped")["count"].sum() / total_count * 100

    # Create custom labels with percentages for both sets of nodes
    clinvar_labels = [f"ClinVar: {label} ({clinvar_percentages[label]:.1f}%)" for label in clinvar_percentages.index]
    am_labels = [f"AM: {label} ({am_percentages[label]:.1f}%)" for label in am_percentages.index]
    node_labels = clinvar_labels + am_labels

    # Create mappings for custom names to indices
    input_mapping = {label: idx for idx, label in enumerate(clinvar_percentages.index)}
    output_mapping = {label: idx + len(input_mapping) for idx, label in enumerate(am_percentages.index)}

    # Update the DataFrame with new mappings
    unival_df["source"] = unival_df["classification_mapped"].map(input_mapping)
    unival_df["target"] = unival_df["pathogenicity_mapped"].map(output_mapping)

    # Create Sankey diagram
    sankey_fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
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

    # Update layout and save Sankey diagram as PNG
    sankey_fig.update_layout(
        title_text=f"{gene_id} ClinVar vs. AlphaMissense Classification",
        font_size=15,
        autosize=False,
        width=1200,
        height=800,
    )
    sankey_path = out_dir / f"{gene_id}_sankey_diagram.png"
    pio.write_image(sankey_fig, file=sankey_path, format="png")
    print(f"Sankey diagram stored at {sankey_path}")
