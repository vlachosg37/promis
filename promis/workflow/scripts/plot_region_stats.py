"""Plot scatter, heatmap, and cytoband instability visualizations per sample.

Inputs
------
- CSV from ``analyze_MSI_distribution.py`` with locus-level instability calls.
- Cytoband annotation file for chromosomal plotting.

Outputs
-------
- Scatter, heatmap, and cytoband PDF plots summarizing regional MSI patterns.
"""

import argparse
import logging
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def plot_cytoband_instability(distr_df, cytobands_file, output_cytoband_file):
    """
    Plot cytoband enrichment analysis based on instability (% Unstable Regions).
    Grouped by chromosome arms (p/q) instead of exact cytobands.
    """
    logger.info(f"Loading cytoband data from {cytobands_file}")
    if distr_df.empty:
        plt.figure(figsize=(8, 4))
        plt.title("Cytoband Arm Instability (% Unstable Regions)")
        plt.text(0.5, 0.5, "No callable regions", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_cytoband_file, dpi=600)
        plt.close()
        return

    cytobands_df = pd.read_csv(
        cytobands_file, sep="\t", names=["Chromosome", "Start", "End", "Band", "Stain"]
    )

    # Ensure 'chr' prefix consistency in cytobands
    cytobands_df["Chromosome"] = cytobands_df["Chromosome"].str.replace("chr", "")

    # Assign p/q chromosome arms to each region
    cytoband_mapping = []

    for _, region in distr_df.iterrows():
        chrom = region["Chromosome"]
        start = region["Region_Start"]
        end = region["Region_End"]

        cytoband_match = cytobands_df[
            (cytobands_df["Chromosome"] == chrom)
            & (cytobands_df["Start"] <= start)
            & (cytobands_df["End"] >= end)
        ]

        if not cytoband_match.empty:
            cytoband_band = cytoband_match.iloc[0]["Band"]
            cytoband_arm = cytoband_band[0]  # 'p' or 'q'
            cytoband_mapping.append(f"{chrom}{cytoband_arm}")
        else:
            cytoband_mapping.append("Unknown")

    distr_df["Cytoband_Arm"] = cytoband_mapping

    # Aggregate cytoband arm statistics
    cytoband_stats = (
        distr_df.groupby("Cytoband_Arm")
        .agg(
            Total_Regions=("MSI_Status", "size"),
            Unstable_Regions=("MSI_Status", lambda x: (x == "Unstable").sum()),
        )
        .reset_index()
    )

    # Normalize instability
    cytoband_stats["% Unstable"] = (
        cytoband_stats["Unstable_Regions"] / cytoband_stats["Total_Regions"]
    ) * 100

    # Remove unknown cytobands if desired
    cytoband_stats = cytoband_stats[cytoband_stats["Cytoband_Arm"] != "Unknown"]

    # Custom sorting function for chromosome arms (e.g., 1p, 1q, 2p, 2q, Xp, Xq)
    def sort_chromosome_arms(arm):
        chrom_part = re.findall(r"(\d+|X|Y)", arm)[0]  # Extract number or X/Y
        arm_part = arm[-1]  # p or q
        chrom_key = int(chrom_part) if chrom_part.isdigit() else {"X": 23, "Y": 24}[chrom_part]
        return (chrom_key, arm_part)

    # Sort cytobands naturally for plotting
    cytoband_stats["Sort_Key"] = cytoband_stats["Cytoband_Arm"].apply(sort_chromosome_arms)
    cytoband_stats = cytoband_stats.sort_values(by="Sort_Key")

    plt.figure(figsize=(12, 6))
    sns.barplot(x="Cytoband_Arm", y="% Unstable", data=cytoband_stats, color="darkred")
    plt.xticks(rotation=90, fontsize=8)
    plt.ylabel("% Unstable Regions")
    plt.xlabel("Chromosome Arm (p/q)")
    plt.title("Cytoband Arm Instability (% Unstable Regions)")
    plt.tight_layout()

    # Save plot
    plt.savefig(output_cytoband_file, dpi=600)
    plt.close()
    logger.info(f"Cytoband arm instability plot saved to {output_cytoband_file}")


def plot_mean_stddev_variation(
    distr_file,
    output_scatter_file,
    output_heatmap_file,
    cytoband_file=None,
    output_cytoband_file=None,
):
    # Load data
    logger.info(f"Loading distribution data from {distr_file}")
    distr_df = pd.read_csv(distr_file)
    logger.info(f"Distribution file columns: {distr_df.columns.tolist()}")
    distr_df = distr_df[distr_df["Chromosome"] != "Summary"].copy()

    if distr_df.empty:
        for output_file, title in (
            (output_scatter_file, "Repeat Length with Variability"),
            (output_heatmap_file, "% Deviating Reads Across Regions"),
        ):
            plt.figure(figsize=(8, 4))
            plt.title(title)
            plt.text(0.5, 0.5, "No callable regions", ha="center", va="center")
            plt.axis("off")
            plt.tight_layout()
            plt.savefig(output_file, dpi=600)
            plt.close()
        if cytoband_file and output_cytoband_file:
            plot_cytoband_instability(distr_df, cytoband_file, output_cytoband_file)
        return

    # Remove 'chr' for consistent chromosome grouping
    distr_df["Chromosome"] = distr_df["Chromosome"].str.replace("chr", "")

    # Handle numeric and non-numeric chromosomes for proper sorting
    def chromosome_sort_key(chrom):
        try:
            return int(chrom)
        except ValueError:
            # Handle 'X', 'Y', and other non-numeric chromosomes
            return float("inf") if chrom == "X" else float("inf") + 1

    distr_df["Chromosome_Sort"] = distr_df["Chromosome"].apply(chromosome_sort_key)
    distr_df = distr_df.sort_values(by=["Chromosome_Sort", "Region_Start"])

    # Create a label combining chromosome and region
    distr_df["Region_Label"] = (
        distr_df["Chromosome"]
        + ":"
        + distr_df["Region_Start"].astype(str)
        + "-"
        + distr_df["Region_End"].astype(str)
    )

    # === Scatter plot with StdDev and % Deviating Reads === #
    plt.figure(figsize=(14, 6))
    scatter = plt.scatter(
        x=range(len(distr_df)),
        y=distr_df["Mean"],
        c=distr_df["% Deviating Reads"],
        cmap="Reds",
        edgecolor="black",
        s=50,
    )

    # Add error bars for StdDev
    plt.errorbar(
        x=range(len(distr_df)),
        y=distr_df["Mean"],
        yerr=distr_df["StdDev"],
        fmt="none",
        ecolor="gray",
        alpha=0.5,
        capsize=3,
    )

    # Add color bar
    plt.colorbar(scatter, label="% Deviating Reads")

    # X-axis labeling, grouped by chromosome
    xtick_labels = []
    prev_chromosome = None
    for _, row in distr_df.iterrows():
        chrom = row["Chromosome"]
        if chrom != prev_chromosome:
            xtick_labels.append(chrom)
            prev_chromosome = chrom
        else:
            xtick_labels.append("")

    plt.xticks(range(len(xtick_labels)), xtick_labels, rotation=90, fontsize=8)
    plt.ylabel("Mean Observed Repeat Length")
    plt.xlabel("Regions (Grouped by Chromosome)")
    plt.title("Repeat Length with Variability (StdDev) and % Deviating Reads (Color)")
    plt.tight_layout()

    # Save scatter plot
    plt.savefig(output_scatter_file, dpi=600)
    plt.close()
    logger.info(f"Scatter plot saved to {output_scatter_file}")

    # === Heatmap for % Deviating Reads === #
    plt.figure(figsize=(14, 8))
    heatmap_data = distr_df.pivot_table(
        index="Chromosome", columns="Region_Start", values="% Deviating Reads", fill_value=0
    )

    chrom_order = sorted(heatmap_data.index, key=chromosome_sort_key)
    heatmap_data = heatmap_data.reindex(chrom_order)

    sns.heatmap(
        heatmap_data,
        cmap="Reds",
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "% Deviating Reads"},
    )

    plt.title("% Deviating Reads Across Regions")
    plt.xlabel("Region Start")
    plt.ylabel("Chromosome")
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save heatmap plot
    plt.savefig(output_heatmap_file, dpi=600)
    plt.close()
    logger.info(f"Heatmap plot saved to {output_heatmap_file}")

    # === Cytoband Plot (Optional) === #
    if cytoband_file and output_cytoband_file:
        plot_cytoband_instability(distr_df, cytoband_file, output_cytoband_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot mean repeat length (scatter) with StdDev error bars and heatmap by % Deviating Reads."
    )
    parser.add_argument("-d", "--distr", required=True, help="Path to the region summary CSV file.")
    parser.add_argument(
        "-s", "--scatter_output", required=True, help="Path to save the scatter plot."
    )
    parser.add_argument(
        "-m", "--heatmap_output", required=True, help="Path to save the heatmap plot."
    )
    parser.add_argument(
        "-c",
        "--cytoband",
        required=False,
        help="Path to cytoband file for cytoband instability plot.",
    )
    parser.add_argument(
        "-cyto_output",
        "--cytoband_output",
        required=False,
        help="Path to save the cytoband instability plot.",
    )

    args = parser.parse_args()
    plot_mean_stddev_variation(
        args.distr,
        args.scatter_output,
        args.heatmap_output,
        cytoband_file=args.cytoband,
        output_cytoband_file=args.cytoband_output,
    )
