"""Classify MSI loci as stable or unstable using deviating read statistics.

Inputs
------
- CSV from ``analyze_MSI_lengths.py`` with per-locus repeat length estimates.

Outputs
-------
- CSV summarizing per-locus instability calls, MSI scores, and MSI status text
  used by downstream plotting and aggregation steps.
"""

import argparse
import logging

import numpy as np
import pandas as pd
from rich.progress import Progress
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


def load_data(input_file):
    """
    Load the input CSV file containing repeat analysis results.

    Args:
        input_file (str): Path to the input CSV file.

    Returns:
        pd.DataFrame: DataFrame with the loaded data.
    """
    logger.info(f"Loading data from: {input_file}")
    return pd.read_csv(input_file)


def fit_gmm_and_get_modeled_length(
    total_lengths, min_fraction=0.25, balance_tolerance=0.025, reg_covar=1e-2, n_init=5
):
    """
    Fit a Gaussian Mixture Model (GMM) to determine the stable allele length(s),
    selecting the best number of components using Bayesian Information Criterion (BIC).

    Args:
        total_lengths (pd.Series): Observed repeat lengths.
        min_fraction (float): Minimum fraction of reads the secondary allele
            must represent.
        balance_tolerance (float): Allowed fractional difference between the
            two major alleles for calling a polymorphism.

    Returns:
        list: Estimated stable repeat lengths.
    """
    MIN_SEPARATION = 1.0
    if len(total_lengths) < 2:
        # Need at least two reads to consider a polymorphism
        logger.warning("Not enough data points for GMM. Using most common length.")
        return [total_lengths.mode().iloc[0]]

    # Determine how many distinct lengths are present. If all reads are the same
    # length, skip GMM fitting entirely.
    unique_lengths = np.unique(total_lengths)
    if len(unique_lengths) == 1:
        logger.debug("Only one unique length detected. Returning the consensus length.")
        return [int(unique_lengths[0])]

    lengths_array = total_lengths.values.reshape(-1, 1)

    best_bic = np.inf
    best_model = None
    best_n_components = 1

    # Limit the number of components to the available distinct lengths and total
    # reads to avoid fitting more clusters than data points.
    max_components = min(3, len(lengths_array), len(unique_lengths))

    # Try different numbers of components and pick the best based on BIC
    for n in range(1, max_components + 1):
        gmm = GaussianMixture(n_components=n, random_state=42, reg_covar=reg_covar, n_init=n_init)
        gmm.fit(lengths_array)
        bic = gmm.bic(lengths_array)

        logger.debug(f"GMM with {n} components: BIC={bic}")

        if bic < best_bic:
            best_bic = bic
            best_model = gmm
            best_n_components = n

    logger.info(f"Selected GMM with {best_n_components} components")

    # Extract means and weights
    gmm_means = best_model.means_.flatten()
    gmm_weights = best_model.weights_.flatten()

    # Sort by frequency
    sorted_indices = np.argsort(gmm_weights)[::-1]
    sorted_means = np.round(gmm_means[sorted_indices]).astype(int)  # Ensure integer repeat lengths
    sorted_weights = gmm_weights[sorted_indices]

    stable_lengths = [sorted_means[0]]  # Always take the most frequent peak

    # Consider a second allele only if it is sufficiently separated, frequent,
    # and has a balanced read fraction relative to the primary allele
    if best_n_components > 1:
        separation = abs(sorted_means[0] - sorted_means[1])
        second_fraction = sorted_weights[1]
        ratio = second_fraction / sorted_weights[0] if sorted_weights[0] > 0 else 0
        logger.debug(
            f"Second component separation={separation}, fraction={second_fraction:.3f}, ratio={ratio:.3f}"
        )
        if (
            separation >= MIN_SEPARATION
            and second_fraction >= min_fraction
            and abs(ratio - 1.0) <= balance_tolerance
        ):
            logger.info(f"Polymorphic sample detected: {sorted_means[0]} & {sorted_means[1]}")
            stable_lengths.append(sorted_means[1])
        else:
            logger.debug("Secondary component ignored due to low support, separation, or imbalance")

    return stable_lengths


def calculate_distribution_stats(
    group,
    call_by="count",
    min_dev_reads=5,
    min_dev_percent=1.0,
    use_GMM=False,
    min_length_percent=0.0,
    balance_tolerance=0.05,
    min_fraction=0.25,
):
    """
    Calculate distribution statistics for a microsatellite region.

    Args:
        group (pd.DataFrame): DataFrame for a specific microsatellite region.
        balance_tolerance (float): Allowed fractional difference between the two
        major alleles when using the GMM approach.
        min_fraction (float): Minimum fraction of reads the secondary allele must
        represent when using the GMM approach.

    Returns:
        dict: A dictionary of computed statistics.
    """
    total_lengths = group["Total_Length_With_Extensions"]

    # Compute basic statistics
    mean_length = total_lengths.mean()
    median_length = total_lengths.median()
    std_length = total_lengths.std()
    min_length = total_lengths.min()
    max_length = total_lengths.max()

    # Detect stable alleles dynamically
    if use_GMM:
        stable_lengths = fit_gmm_and_get_modeled_length(
            total_lengths,
            min_fraction=min_fraction,
            balance_tolerance=balance_tolerance,
        )
        if len(stable_lengths) == 1:
            # Fall back to the median when only one GMM peak is detected
            stable_lengths = [round(total_lengths.median())]
    else:
        stable_lengths = [round(total_lengths.median())]

    # Define deviation threshold in base pairs. Reads differing by at least one base from all stable lengths are considered deviating.
    deviation_threshold = 1

    # Identify deviating reads (outside threshold)
    deviating_reads = total_lengths[
        ~total_lengths.isin(stable_lengths)
        & (
            (total_lengths <= min(stable_lengths) - deviation_threshold)
            | (total_lengths >= max(stable_lengths) + deviation_threshold)
        )
    ]

    if min_length_percent > 0:
        # Determine the maximum count among stable lengths as reference
        length_counts = total_lengths.value_counts()
        stable_counts = length_counts[length_counts.index.isin(stable_lengths)]
        reference_count = stable_counts.max() if not stable_counts.empty else len(total_lengths)
        cutoff = reference_count * (min_length_percent / 100.0)

        # Filter out low frequency lengths from deviating reads
        dev_counts = deviating_reads.value_counts()
        valid_lengths = dev_counts[dev_counts >= cutoff].index
        deviating_reads = deviating_reads[deviating_reads.isin(valid_lengths)]

    # Calculate percentage of deviating reads
    percent_deviating = (len(deviating_reads) / len(total_lengths)) * 100

    # Determine MSI status based on user-defined strategy
    if call_by == "count":
        msi_status = "Unstable" if len(deviating_reads) >= min_dev_reads else "Stable"
    elif call_by == "percent":
        msi_status = "Unstable" if percent_deviating >= min_dev_percent else "Stable"
    elif call_by == "both":
        msi_status = (
            "Unstable"
            if (len(deviating_reads) >= min_dev_reads and percent_deviating >= min_dev_percent)
            else "Stable"
        )
    else:
        raise ValueError("Invalid value for --call_by. Choose from: 'count', 'percent', or 'both'.")

    return {
        "Mean": mean_length,
        "Median": median_length,
        "Stable_Alleles": stable_lengths,
        "StdDev": std_length,
        "Min": min_length,
        "Max": max_length,
        "% Deviating Reads": percent_deviating,
        "Deviating_Reads": len(deviating_reads),
        "MSI_Status": msi_status,
        "Expected_Length_Reads": len(total_lengths),
    }


def analyze_distribution(
    input_file,
    output_path,
    use_mode=False,
    use_GMM=False,
    call_by="count",
    min_dev_reads=5,
    min_dev_percent=1.0,
    min_length_percent=0.0,
    balance_tolerance=0.05,
    min_fraction=0.25,
    min_total_reads=50,
):
    """
    Analyze the distribution of repeat lengths within each microsatellite region.

    Args:
        input_file (str): Path to the input CSV file.
        output_path (str): Path to save the output CSV file with statistics.
        use_mode (bool): If True, uses mode instead of median for MSI classification.
        balance_tolerance (float): Allowed fractional difference between primary
        and secondary alleles when using the GMM approach.
        min_fraction (float): Minimum fraction of reads the secondary allele must
        represent when using the GMM approach.
        min_total_reads (int): Skip regions with fewer total reads than this threshold.
    """
    # Load input data
    data = load_data(input_file)

    # Filter for reads where Context_Match is 'Pass'
    data = data[data["Context_Match"] == "Pass"]

    # Ensure Total_Length_With_Extensions is numeric and drop invalid rows
    data["Total_Length_With_Extensions"] = pd.to_numeric(
        data["Total_Length_With_Extensions"], errors="coerce"
    )
    data = data.dropna(subset=["Total_Length_With_Extensions"])

    # Sort data by Chromosome, Region_Start, and Region_End
    data["Chromosome_Sort"] = data["Chromosome"].str.extract(r"(\d+)").fillna(0).astype(int)
    data = data.sort_values(by=["Chromosome_Sort", "Region_Start", "Region_End"])
    data["Chromosome"] = data["Chromosome_Sort"].apply(lambda x: f"chr{x}" if x != 0 else "chrX")

    # Group by Chromosome, Region_Start, and Region_End
    grouped = data.groupby(["Chromosome", "Region_Start", "Region_End"])

    # Collect statistics for each group
    results = []
    with Progress() as progress:
        task = progress.add_task("Analyzing regions", total=grouped.ngroups)
        for (chromosome, start, end), group in grouped:
            if len(group) < min_total_reads:
                logger.info(
                    f"Skipping {chromosome}:{start}-{end} due to low coverage ({len(group)} reads)"
                )
                progress.advance(task)
                continue
            stats = calculate_distribution_stats(
                group,
                call_by=call_by,
                min_dev_reads=min_dev_reads,
                min_dev_percent=min_dev_percent,
                use_GMM=use_GMM,
                min_length_percent=min_length_percent,
                balance_tolerance=balance_tolerance,
                min_fraction=min_fraction,
            )

            stats.update(
                {
                    "Chromosome": chromosome,
                    "Region_Start": start,
                    "Region_End": end,
                    "Expected_Length": group["Expected_Length"].iloc[0],
                }
            )
            results.append(stats)
            progress.advance(task)

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # If no regions met the coverage threshold, create an empty results file so
    # downstream steps can proceed.
    if results_df.empty:
        logger.info("Sample skipped: no regions passed the coverage threshold.")
        empty_cols = [
            "Chromosome",
            "Region_Start",
            "Region_End",
            "Expected_Length",
            "Expected_Length_Reads",
            "Deviating_Reads",
            "Mean",
            "Median",
            "Stable_Alleles",
            "StdDev",
            "Min",
            "Max",
            "% Deviating Reads",
            "MSI_Status",
        ]
        pd.DataFrame(columns=empty_cols).to_csv(output_path, index=False)
        return

    # Reorder columns to place Chromosome, Region_Start, Region_End, and Expected_Length first
    cols = [
        "Chromosome",
        "Region_Start",
        "Region_End",
        "Expected_Length",
        "Expected_Length_Reads",
        "Deviating_Reads",
    ] + [
        col
        for col in results_df.columns
        if col
        not in [
            "Chromosome",
            "Region_Start",
            "Region_End",
            "Expected_Length",
            "Expected_Length_Reads",
            "Deviating_Reads",
        ]
    ]
    results_df = results_df[cols]

    # Calculate overall summary stats for the sample
    total_regions = len(results_df)
    unstable_regions = (results_df["MSI_Status"] == "Unstable").sum()
    percent_unstable = (unstable_regions / total_regions) * 100 if total_regions > 0 else 0

    # Add summary row
    summary_row = {
        "Chromosome": "Summary",
        "Region_Start": None,
        "Region_End": None,
        "Expected_Length": None,
        "Expected_Length_Reads": None,
        "Deviating_Reads": None,
        "Mean": None,
        "Median": None,
        "StdDev": None,
        "Min": None,
        "Max": None,
        "% Deviating Reads": None,
        "MSI_Status": f"{percent_unstable:.2f}% Unstable ({unstable_regions}/{total_regions} regions)",
    }

    results_df = pd.concat([results_df, pd.DataFrame([summary_row])], ignore_index=True)

    # Ensure Chromosome is sorted naturally: chr1, chr2, ..., chrX, chrY
    def chromosome_sort_key(chrom):
        try:
            return int(chrom.replace("chr", ""))
        except ValueError:
            return float("inf") if chrom == "chrX" else float("inf") + 1

    # Apply sorting before saving
    results_df["Chromosome"] = results_df["Chromosome"].astype(str)  # Ensure string format
    results_df = results_df.sort_values(by=["Chromosome"], key=lambda x: x.map(chromosome_sort_key))

    # Save results to CSV
    results_df.to_csv(output_path, index=False)

    logger.info(f"Analysis results saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze intra-sample repeat length distributions for MSI detection."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output CSV file.")
    parser.add_argument(
        "--use_mode", action="store_true", help="Use mode instead of median for MSI classification."
    )
    parser.add_argument(
        "--use_GMM",
        action="store_true",
        help="Use Gaussian Mixture Model (GMM) for MSI classification.",
    )
    parser.add_argument(
        "--call_by",
        choices=["count", "percent", "both"],
        default="count",
        help="Method to determine MSI status: 'count', 'percent', or 'both' (default: count)",
    )
    parser.add_argument(
        "--min_dev_reads",
        type=int,
        default=5,
        help="Minimum number of deviating reads to call a region unstable (default: 5)",
    )
    parser.add_argument(
        "--min_dev_percent",
        type=float,
        default=1.0,
        help="Minimum percent of deviating reads to call a region unstable (default: 1.0)",
    )
    parser.add_argument(
        "--min_length_percent",
        type=float,
        default=0.0,
        help="Minimum percent of the main stable length that a deviating length must represent to be considered (default: 0.0)",
    )
    parser.add_argument(
        "--balance_tolerance",
        type=float,
        default=0.05,
        help="Allowed fractional difference between primary and secondary alleles when using GMM (default: 0.05)",
    )
    parser.add_argument(
        "--min_fraction",
        type=float,
        default=0.25,
        help="Minimum fraction of reads the secondary allele must represent when using GMM (default: 0.25)",
    )
    parser.add_argument(
        "--min_total_reads",
        type=int,
        default=50,
        help="Exclude regions with total reads below this number from scoring (default: 50)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    parser.add_argument("--info", action="store_true", help="Enable info-level logging.")

    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO if args.info else logging.WARNING
    logging.basicConfig(level=log_level, format="%(levelname)s:%(message)s")
    logger.setLevel(log_level)

    analyze_distribution(
        args.input,
        args.output,
        use_mode=args.use_mode,
        use_GMM=args.use_GMM,
        call_by=args.call_by,
        min_dev_reads=args.min_dev_reads,
        min_dev_percent=args.min_dev_percent,
        min_length_percent=args.min_length_percent,
        balance_tolerance=args.balance_tolerance,
        min_fraction=args.min_fraction,
        min_total_reads=args.min_total_reads,
    )


if __name__ == "__main__":
    main()
