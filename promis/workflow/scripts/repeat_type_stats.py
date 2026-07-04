"""Summarize instability patterns by repeat unit length and type.

Inputs
------
- CSV from ``analyze_MSI_lengths.py`` with per-locus repeat length estimates.

Outputs
-------
- Summary CSV aggregating instability metrics by repeat unit length.
- Optional barplot and scatterplot PDFs highlighting unstable repeat classes.
"""

import argparse
import logging
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_repeat_length(repeat):
    matches = re.findall(r"\((.*?)\)(\d+)", str(repeat))
    if matches:
        return sum(int(length) for _, length in matches)
    return None


def classify_repeat_type(repeat):
    matches = re.findall(r"\((.*?)\)\d+", str(repeat))
    unique_units = set(matches)
    if len(unique_units) == 1:
        unit_length = len(next(iter(unique_units)))
        if unit_length == 1:
            return "Mononucleotide"
        elif unit_length == 2:
            return "Dinucleotide"
        elif unit_length == 3:
            return "Trinucleotide"
    return "Complex"


def analyze_repeat_types(
    lengths_file,
    output_dir,
    plot_frequency,
    plot_instability,
    plot_composition,
    plot_length_vs_instability,
    sample,
):
    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"Loading lengths data from {lengths_file}")
    lengths_df = pd.read_csv(lengths_file)
    lengths_df = lengths_df[lengths_df["Context_Match"] == "Pass"]

    lengths_df["Repeat_Length"] = lengths_df["Expected_Repeat"].apply(parse_repeat_length)
    lengths_df["Repeat_Type"] = lengths_df["Expected_Repeat"].apply(classify_repeat_type)
    lengths_df["Repeat_Unit"] = lengths_df["Expected_Repeat"].apply(
        lambda x: "".join(re.findall(r"\((.*?)\)", str(x)))
    )

    repeat_counts = lengths_df["Expected_Repeat"].value_counts().reset_index()
    repeat_counts.columns = ["Expected_Repeat", "Count"]

    stability_group = (
        lengths_df.groupby("Repeat_Unit")["MSI_Status"]
        .value_counts(normalize=True)
        .unstack(fill_value=0)
    )
    if "Unstable" in stability_group.columns:
        stability_group["% Unstable"] = stability_group["Unstable"] * 100
    else:
        stability_group["% Unstable"] = 0

    summary_df = lengths_df.groupby("Repeat_Unit").agg(
        Count=("Repeat_Unit", "count"),
        Min_Length=("Repeat_Length", "min"),
        Max_Length=("Repeat_Length", "max"),
        Median_Length=("Repeat_Length", "median"),
        Std_Length=("Repeat_Length", "std"),
    )

    summary_df = summary_df.merge(
        stability_group["% Unstable"], left_index=True, right_index=True, how="left"
    ).fillna(0)
    summary_df.reset_index(inplace=True)
    summary_df.to_csv(f"{output_dir}/{sample}_repeat_type_summary.csv", index=False)

    if plot_frequency:
        plt.figure(figsize=(10, 6))
        sns.barplot(x="Count", y="Expected_Repeat", data=repeat_counts.head(20))
        plt.title("Top 20 Repeat Types by Frequency")
        plt.xlabel("Count")
        plt.ylabel("Repeat Type")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{sample}_repeat_type_frequency_barplot.pdf")
        plt.close()

    if plot_instability:
        stability_sorted = (
            stability_group["% Unstable"].sort_values(ascending=False).head(20).reset_index()
        )
        plt.figure(figsize=(10, 6))
        sns.barplot(x="% Unstable", y="Repeat_Unit", data=stability_sorted)
        plt.title("Top 20 Repeat Units by % Unstable")
        plt.xlabel("% Unstable")
        plt.ylabel("Repeat Unit")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{sample}_repeat_unit_instability_barplot.pdf")
        plt.close()

    if plot_composition:
        type_counts = lengths_df["Repeat_Type"].value_counts(normalize=True) * 100
        plt.figure(figsize=(6, 6))
        type_counts.plot(
            kind="pie", autopct="%1.1f%%", startangle=90, colors=sns.color_palette("pastel")
        )
        plt.title("Repeat Type Composition")
        plt.ylabel("")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{sample}_repeat_type_composition_piechart.pdf")
        plt.close()

        plt.figure(figsize=(8, 6))
        sns.barplot(x=type_counts.index, y=type_counts.values, palette="pastel")
        plt.title("Repeat Type Composition")
        plt.xlabel("Repeat Type")
        plt.ylabel("Percentage")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{sample}_repeat_type_composition_barplot.pdf")
        plt.close()

    if plot_length_vs_instability:
        length_group = (
            lengths_df.groupby("Repeat_Length")["MSI_Status"]
            .value_counts(normalize=True)
            .unstack(fill_value=0)
        )
        if "Unstable" in length_group.columns:
            length_group["% Unstable"] = length_group["Unstable"] * 100
        else:
            length_group["% Unstable"] = 0

        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=length_group.index, y=length_group["% Unstable"])
        plt.title("Repeat Length vs. % Unstable")
        plt.xlabel("Repeat Length")
        plt.ylabel("% Unstable")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{sample}_repeat_length_vs_instability_scatterplot.pdf")
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze repeat types and plot statistics.")
    parser.add_argument("-l", "--lengths", required=True, help="Path to the lengths CSV file")
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to save plots and summary"
    )
    parser.add_argument(
        "--plot_frequency", action="store_true", help="Generate repeat type frequency plot"
    )
    parser.add_argument(
        "--plot_instability", action="store_true", help="Generate % unstable by repeat type plot"
    )
    parser.add_argument(
        "--plot_composition",
        action="store_true",
        help="Generate repeat type composition pie and bar plots",
    )
    parser.add_argument(
        "--plot_length_vs_instability",
        action="store_true",
        help="Generate repeat length vs. % unstable plot",
    )
    parser.add_argument("--sample", required=True, help="Sample name for output files")
    args = parser.parse_args()

    analyze_repeat_types(
        args.lengths,
        args.output_dir,
        args.plot_frequency,
        args.plot_instability,
        args.plot_composition,
        args.plot_length_vs_instability,
        args.sample,
    )
