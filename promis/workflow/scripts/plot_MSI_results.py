"""Create per-sample MSI status barplots from distribution summaries.

Inputs
------
- CSV from ``analyze_MSI_distribution.py`` containing MSI_Status and per-locus
  calls for a single sample.

Outputs
-------
- PDF barplot showing counts of stable vs. unstable loci for the sample.
"""

import argparse
import logging

import matplotlib.pyplot as plt
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def load_analysis_results(input_file):
    """
    Load the analysis results CSV file.

    Args:
        input_file (str): Path to the analysis results file.

    Returns:
        pd.DataFrame: DataFrame with analysis results.
    """
    logger.info(f"Loading analysis results from: {input_file}")
    return pd.read_csv(input_file)


def calculate_msi_percentages(results_df):
    """
    Calculate the percentage of stable and unstable regions.

    Args:
        results_df (pd.DataFrame): DataFrame with analysis results.

    Returns:
        tuple: Percentages of unstable and stable regions.
    """
    total_regions = len(results_df)
    unstable_count = len(results_df[results_df["MSI_Status"] == "Unstable"])
    stable_count = total_regions - unstable_count

    percent_unstable = (unstable_count / total_regions) * 100
    percent_stable = (stable_count / total_regions) * 100

    return percent_unstable, percent_stable


def plot_msi_status(percent_unstable, percent_stable, output_file, sample_name):
    """
    Create a stacked barplot showing the percentage of unstable and stable regions.

    Args:
        percent_unstable (float): Percentage of unstable regions.
        percent_stable (float): Percentage of stable regions.
        output_file (str): Path to save the barplot.
        sample_name (str): Name of the sample to include in the title.
    """
    logger.info("Creating stacked barplot for MSI status percentages.")

    categories = [""]
    unstable = [percent_unstable]
    stable = [percent_stable]

    fig, ax = plt.subplots(figsize=(6, 5))

    ax.bar(categories, unstable, label="Unstable", color="darkred")
    ax.bar(categories, stable, bottom=unstable, label="Stable", color="lightgrey")

    # Annotate only the unstable percentage inside the red bar
    ax.text(
        0,
        unstable[0] / 2,
        f"{percent_unstable:.1f}%",
        ha="center",
        va="center",
        color="white",
        fontsize=10,
    )

    # Formatting improvements
    ax.set_ylabel("Percentage of Regions", labelpad=10)
    ax.set_title(f"MSI Status: {sample_name}")
    ax.set_ylim(0, 100)
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(output_file, dpi=600)
    logger.info(f"Barplot saved to: {output_file}")
    plt.close()


def main(input_file, output_file):
    """
    Main function to load data, calculate percentages, and create the barplot.

    Args:
        input_file (str): Path to the input analysis results file.
        output_file (str): Path to save the barplot.
    """
    # Load data
    results_df = load_analysis_results(input_file)

    # Extract sample name from the input file path
    sample_name = input_file.split("/")[-1].split("_marked")[0]

    # Calculate percentages
    percent_unstable, percent_stable = calculate_msi_percentages(results_df)

    # Create barplot
    plot_msi_status(percent_unstable, percent_stable, output_file, sample_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a barplot for MSI status distribution.")
    parser.add_argument(
        "-i", "--input", required=True, help="Path to the input CSV file with analysis results."
    )
    parser.add_argument("-o", "--output", required=True, help="Path to save the barplot.")

    args = parser.parse_args()

    main(args.input, args.output)
