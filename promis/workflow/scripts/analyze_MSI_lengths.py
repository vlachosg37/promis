"""Quantify repeat lengths for microsatellite loci using extracted reads.

Inputs
------
- CSV produced by ``extract_MSI_sequences.py`` with read-level alignments per
  locus.
- Reference MSI loci metadata CSV.

Outputs
-------
- Per-sample CSV reporting observed repeat lengths, consensus length estimates,
  and supporting read counts per locus.
"""

import argparse
import logging
import re
from difflib import SequenceMatcher

import numpy as np
import pandas as pd
from numba import njit
from numba.typed import List
from rich.progress import Progress

logger = logging.getLogger(__name__)


def load_extracted_reads(extracted_csv):
    logger.info(f"Loading extracted reads from: {extracted_csv}")
    reads_df = pd.read_csv(extracted_csv)
    return reads_df


def load_reference_data(reference_csv):
    logger.info(f"Loading reference context data from: {reference_csv}")
    reference_df = pd.read_csv(reference_csv)
    return reference_df[["Chromosome", "Start", "End", "Upstream_Context", "Downstream_Context"]]


def parse_expected_repeat(expected_repeat):
    """
    Parse repeat sequences like (T)9 or (T)7(C)7 into a list of (unit, count) tuples.
    """
    if not isinstance(expected_repeat, str) or pd.isna(expected_repeat):
        logger.warning(f"Invalid Expected_Repeat value: {expected_repeat}")
        return []

    pattern = re.findall(r"\((\w+)\)(\d+)", expected_repeat)
    if not pattern:
        logger.warning(f"Unexpected repeat pattern: {expected_repeat}")
        return []

    return [(unit, int(count)) for unit, count in pattern]


def find_approximate_match(sequence, subsequence, min_ratio=0.8):
    seq_len = len(sequence)
    sub_len = len(subsequence)
    for i in range(seq_len - sub_len + 1):
        window = sequence[i : i + sub_len]
        match_ratio = SequenceMatcher(None, window, subsequence).ratio()
        if match_ratio >= min_ratio:
            return i, i + sub_len
    return -1, -1


def calculate_repeat_start_position(row):
    try:
        repeat_start = int(row["Repeat_Coordinates"].split(":")[1].split("-")[0])
        read_start = int(row["Read_Start"])
        position_in_read = repeat_start - read_start
        if 0 <= position_in_read < len(row["Read_Sequence"]):
            logger.info(f"Repeat starts at position {position_in_read} in read: {row['Read_Name']}")
            return position_in_read
        else:
            logger.warning(
                f"Repeat start {repeat_start} not within read bounds ({row['Read_Start']}-{row['Read_End']})"
            )
            return -1
    except Exception as e:
        logger.warning(f"Failed to calculate repeat start position for row {row['Read_Name']}: {e}")
        return -1


def find_repeat_in_sequence(sequence, repeat_sequence, approximate=False):
    if approximate:
        return find_approximate_match(sequence, repeat_sequence)
    start_index = sequence.find(repeat_sequence)
    end_index = start_index + len(repeat_sequence) if start_index != -1 else -1
    return (start_index, end_index)


@njit
def find_repeat_run(sequence, repeat_unit, min_repeats=3):
    """
    Locate the longest run of a repeat unit in a sequence.
    Args:
        sequence (np.ndarray): Array representation of the sequence.
        repeat_unit (np.ndarray): Array of the repeat unit.
        min_repeats (int): Minimum number of repeats to consider a match.

    Returns:
        tuple: ((start, end), number_of_repeats).
    """
    unit_len = repeat_unit.shape[0]
    max_run_start = -1
    max_run_end = -1
    max_run_repeats = 0
    seq_len = sequence.shape[0]

    i = 0
    while i <= seq_len - unit_len:
        match = True
        for j in range(unit_len):
            if sequence[i + j] != repeat_unit[j]:
                match = False
                break
        if match:
            run_start = i
            run_repeats = 0

            while i <= seq_len - unit_len:
                inner_match = True
                for j in range(unit_len):
                    if sequence[i + j] != repeat_unit[j]:
                        inner_match = False
                        break
                if not inner_match:
                    break
                run_repeats += 1
                i += unit_len

            run_end = run_start + run_repeats * unit_len

            if run_repeats >= min_repeats and run_repeats > max_run_repeats:
                max_run_start = run_start
                max_run_end = run_end
                max_run_repeats = run_repeats
        else:
            i += 1

    return (max_run_start, max_run_end), max_run_repeats


@njit
def find_complex_repeat_run(sequence, repeat_components, min_fraction=0.67):
    """
    Locate the best match for a complex repeat (e.g., (T)7(C)7).

    Args:
        sequence (np.ndarray): Sequence to search.
        repeat_components (List): List of (repeat_unit_array, count) tuples.
        min_fraction (float): Minimum fraction of the expected repeat length to consider a valid match.

    Returns:
        tuple: ((start, end), observed_length), or ((-1, -1), 0) if not found.
    """
    expected_repeats = 0
    for _, count in repeat_components:
        expected_repeats += count
    min_repeats_required = int(expected_repeats * min_fraction)

    best_start = -1
    best_end = -1
    best_length = 0
    seq_len = sequence.shape[0]

    i = 0
    while i < seq_len:
        start = i
        observed_length = 0
        match = True

        for unit, expected_count in repeat_components:
            unit_len = unit.shape[0]
            count = 0

            while i <= seq_len - unit_len:
                unit_match = True
                for j in range(unit_len):
                    if sequence[i + j] != unit[j]:
                        unit_match = False
                        break
                if not unit_match:
                    break
                count += 1
                i += unit_len

            observed_length += count * unit_len

            if count < int(expected_count * min_fraction):
                match = False
                break

        if match and observed_length >= min_repeats_required:
            end = i
            if observed_length > best_length:
                best_start = start
                best_end = end
                best_length = observed_length

        if not match:
            i = start + 1

    return (best_start, best_end), best_length


@njit
def check_repeat_extension(sequence, repeat_unit, start, end):
    upstream_count = 0
    downstream_count = 0
    unit_len = repeat_unit.shape[0]
    seq_len = sequence.shape[0]

    i = start - unit_len
    while i >= 0:
        match = True
        for j in range(unit_len):
            if sequence[i + j] != repeat_unit[j]:
                match = False
                break
        if not match:
            break
        upstream_count += 1
        i -= unit_len

    i = end
    while i + unit_len <= seq_len:
        match = True
        for j in range(unit_len):
            if sequence[i + j] != repeat_unit[j]:
                match = False
                break
        if not match:
            break
        downstream_count += 1
        i += unit_len

    total_length = (end - start) + (upstream_count + downstream_count) * unit_len
    return total_length, upstream_count * unit_len, downstream_count * unit_len


def verify_flanking_context(row, repeat_position, observed_length, repeat_units):
    """
    Verify that the upstream and downstream context flanking the repeat matches the expected context.
    Handles both simple and complex repeats.
    """
    try:
        start, end = repeat_position
        if start == -1 or end == -1:
            return "Not Checked"

        upstream_context_read = row["Read_Sequence"][max(0, start - 4) : start]

        # For both simple and complex, adjust downstream based on observed length
        downstream_start = start + observed_length if observed_length != "Not Found" else end

        downstream_context_read = row["Read_Sequence"][downstream_start : downstream_start + 4]

        expected_upstream_context = row.get("Upstream_Context", "")
        expected_downstream_context = row.get("Downstream_Context", "")

        upstream_match = upstream_context_read == expected_upstream_context
        downstream_match = downstream_context_read == expected_downstream_context

        if upstream_match and downstream_match:
            return "Pass"
        return "Fail"
    except Exception as e:
        logger.warning(f"Failed context verification for {row['Read_Name']}: {e}")
        return "Error"


def analyze_repeats(extracted_reads, reference_df, output_file, threshold=0.1):
    # Merge reference context data into extracted reads based on Chromosome, Start, End
    extracted_reads["Chromosome"] = extracted_reads["Chromosome"].astype(str)
    reference_df["Chromosome"] = reference_df["Chromosome"].astype(str)

    extracted_reads = extracted_reads.merge(
        reference_df,
        left_on=["Chromosome", "Region_Start", "Region_End"],
        right_on=["Chromosome", "Start", "End"],
        how="left",
    )

    results = []
    with Progress() as progress:
        task = progress.add_task("Analyzing reads", total=len(extracted_reads))
        for _, row in extracted_reads.iterrows():
            chrom, region_start, region_end, sequence = (
                row["Chromosome"],
                row["Region_Start"],
                row["Region_End"],
                row["Read_Sequence"],
            )
            logger.info(f"Processing region: {chrom}:{region_start}-{region_end}")

            sequence_array = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

            # Parse the repeat components, e.g., [(T, 7), (C, 7)]
            repeat_components_raw = parse_expected_repeat(row["Expected_Repeat"])
            if not repeat_components_raw:
                logger.warning(
                    f"Failed to parse repeat components from Expected_Repeat for region: {chrom}:{region_start}-{region_end}"
                )
                progress.advance(task)
                continue

            # Determine expected repeat length and minimum acceptable repeats (2/3 of expected)
            expected_length = sum(len(unit) * count for unit, count in repeat_components_raw)

            repeat_components = List()
            for unit, count in repeat_components_raw:
                repeat_components.append(
                    (np.frombuffer(unit.encode("ascii"), dtype=np.uint8), count)
                )

            # Locate the best match for the entire complex repeat run
            repeat_position_run, observed_length = find_complex_repeat_run(
                sequence_array, repeat_components, min_fraction=2 / 3
            )

            # Approximate fallback removed for complex repeats – focus on exact run detection
            repeat_position_approx = (-1, -1)

            # Finalize the best position and observed length
            if repeat_position_run != (-1, -1):
                repeat_position = repeat_position_run
            else:
                repeat_position = (-1, -1)
                observed_length = "Not Found"

            if repeat_position == (-1, -1):
                logger.info(f"Repeat not found for region: {chrom}:{region_start}-{region_end}")
                context_match = "Not Checked"
                total_length = "Not Found"
                extensions = {"upstream": "Not Found", "downstream": "Not Found"}
                deviation = "Not Found"
                msi_status = "Not Found in Read"
            else:
                if len(repeat_components_raw) > 1:
                    # Complex repeat - check extensions with both first and last units
                    first_repeat_unit = np.frombuffer(
                        repeat_components_raw[0][0].encode("ascii"), dtype=np.uint8
                    )
                    last_repeat_unit = np.frombuffer(
                        repeat_components_raw[-1][0].encode("ascii"), dtype=np.uint8
                    )

                    tl_first, up_first, down_first = check_repeat_extension(
                        sequence_array, first_repeat_unit, *repeat_position
                    )
                    tl_last, up_last, down_last = check_repeat_extension(
                        sequence_array, last_repeat_unit, *repeat_position
                    )

                    if tl_first > tl_last:
                        total_length, extensions = tl_first, {
                            "upstream": up_first,
                            "downstream": down_first,
                        }
                    else:
                        total_length, extensions = tl_last, {
                            "upstream": up_last,
                            "downstream": down_last,
                        }
                else:
                    # Simple repeat - standard extension check
                    first_repeat_unit = np.frombuffer(
                        repeat_components_raw[0][0].encode("ascii"), dtype=np.uint8
                    )
                    tl, up_ext, down_ext = check_repeat_extension(
                        sequence_array, first_repeat_unit, *repeat_position
                    )
                    total_length = tl
                    extensions = {"upstream": up_ext, "downstream": down_ext}

                deviation = abs(total_length - expected_length)
                msi_status = (
                    "Unstable"
                    if total_length < expected_length or deviation / expected_length >= threshold
                    else "Stable"
                )
                context_match = verify_flanking_context(
                    row, repeat_position, observed_length, repeat_components_raw
                )

            results.append(
                {
                    **row.to_dict(),
                    "Repeat_Position_Exact": repeat_position_run,
                    "Repeat_Position_Approx": repeat_position_approx,
                    "Observed_Length": observed_length,
                    "Total_Length_With_Extensions": total_length,
                    "Extensions": extensions,
                    "Expected_Length": expected_length,
                    "Deviation": deviation,
                    "MSI_Status": msi_status,
                    "Context_Match": context_match,
                }
            )
            progress.advance(task)

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    logger.info(f"Analysis results saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze repeat lengths and determine MSI status.")
    parser.add_argument("-e", "--extracted", required=True, help="Path to the extracted reads CSV.")
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Path to the reference CSV with upstream and downstream contexts.",
    )
    parser.add_argument("-o", "--output", required=True, help="Path to save the analysis results.")
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.1,
        help="Threshold for MSI calling (default: 0.1).",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    parser.add_argument("--info", action="store_true", help="Enable info-level logging.")

    args = parser.parse_args()
    log_level = logging.DEBUG if args.verbose else logging.INFO if args.info else logging.WARNING
    logging.basicConfig(level=log_level, format="%(levelname)s:%(message)s")
    logger.setLevel(log_level)

    extracted_reads = load_extracted_reads(args.extracted)
    reference_data = load_reference_data(args.reference)

    analyze_repeats(extracted_reads, reference_data, args.output, args.threshold)
