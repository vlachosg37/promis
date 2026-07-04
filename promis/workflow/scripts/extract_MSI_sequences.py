"""
===============================================================================
extract_MSI_sequences.py
===============================================================================

Description:
    This script extracts high-quality reads from an alignment file (BAM/CRAM) that fully span
    predefined microsatellite (MSI) repeat regions. It supports consensus BAMs
    with UMIs and enables filtering by base quality, mapping quality, and
    sequence ambiguity.

    If Unique Molecular Identifiers (UMIs) are present in the RX, UR, or UB
    SAM tags, the script will:
      - Extract the UMI per read (recorded as "NA" if absent)
      - Group reads by (Chromosome, Repeat_Coordinates, UMI, Read_Start, Read_End)
      - Collapse each group to a single consensus allele by retaining the
        highest-confidence read supporting the majority sequence
      - Retain the highest-confidence read based on:
            1. Mapping quality
            2. Mean base quality
            3. Read length

Features:
    • Supports configurable base quality and mapping quality thresholds
    • Optional filtering of reads containing ambiguous 'N' bases
    • Automatically resolves chromosome naming mismatches (e.g., 'chr1' vs '1')
    • Logs UMI tag usage and deduplication status
    • Outputs a deduplicated CSV table with MSI region read summaries

Intended Use:
    Designed for pre-processing consensus alignment files in MSI detection pipelines
    that require accurate, per-molecule analysis of microsatellite repeat lengths.

===============================================================================
"""

import argparse
import logging
import os
import statistics

import pandas as pd
import pysam
from rich.progress import Progress, track

# Configure logging
logger = logging.getLogger(__name__)


def load_repeat_coordinates(repeats_file):
    logger.info(f"Loading repeat coordinates from: {repeats_file}")
    repeats_df = pd.read_csv(repeats_file)
    logger.info(f"Columns in repeats_df: {repeats_df.columns}")
    # Sort by chromosome and start to allow sequential BAM access
    repeats_df = repeats_df.sort_values(by=["Chromosome", "Start"]).reset_index(drop=True)
    return repeats_df


def _open_alignment(filepath, reference_path=None, threads=None):
    lower_path = filepath.lower()
    kwargs = {"threads": threads} if threads is not None else {}
    if lower_path.endswith(".cram"):
        if not reference_path:
            raise ValueError("A reference genome path is required when reading CRAM input.")
        return pysam.AlignmentFile(filepath, "rc", reference_filename=reference_path, **kwargs)
    if lower_path.endswith(".bam"):
        return pysam.AlignmentFile(filepath, "rb", **kwargs)
    raise ValueError("Unsupported alignment format. Expected a .bam or .cram file.")


def check_chr_format(alignment_file, repeats_df, reference_path=None):
    repeats_df["Chromosome"] = repeats_df["Chromosome"].astype(str)
    with _open_alignment(alignment_file, reference_path=reference_path) as alignment:
        alignment_chromosomes = alignment.references
    bam_uses_chr = alignment_chromosomes[0].startswith("chr")
    repeats_uses_chr = repeats_df["Chromosome"].iloc[0].startswith("chr")
    if bam_uses_chr and not repeats_uses_chr:
        repeats_df["Chromosome"] = repeats_df["Chromosome"].apply(
            lambda x: f"chr{x}" if not x.startswith("chr") else x
        )
    elif not bam_uses_chr and repeats_uses_chr:
        repeats_df["Chromosome"] = repeats_df["Chromosome"].str.replace("^chr", "", regex=True)
    return repeats_df


def calculate_base_quality_stats(qualities):
    if not qualities:
        return {"Mean_Quality": None, "Median_Quality": None}
    mean_quality = round(statistics.mean(qualities), 2)
    median_quality = round(statistics.median(qualities), 2)
    return {"Mean_Quality": mean_quality, "Median_Quality": median_quality}


def collapse_umi_families(df):
    """Collapse reads sharing UMI and mapping position to a consensus.

    Families are defined by UMI together with the aligned start and end
    coordinates. Within each family the read sequence occurring most often is
    taken as the true allele; ties are broken by mapping quality, mean base
    quality and read length.
    """
    consensus_reads = []
    for _, family in df.groupby(["UMI", "Read_Start", "Read_End"]):
        # Identify the majority allele (most common sequence)
        majority_seq = family["Read_Sequence"].value_counts().idxmax()
        majority_reads = family[family["Read_Sequence"] == majority_seq]
        # Select the highest-confidence read supporting the majority allele
        best_read = majority_reads.sort_values(
            by=["Mapping_Quality", "Mean_Quality", "Read_Length"],
            ascending=False,
        ).iloc[0]
        consensus_reads.append(best_read)
    return pd.DataFrame(consensus_reads)


def extract_reads_from_alignment(
    alignment_file,
    repeats_df,
    output_file,
    bq_threshold,
    mq_threshold,
    keep_n,
    min_reads,
    reference_path=None,
):
    try:
        sample_name = os.path.splitext(os.path.basename(alignment_file))[0]
        logger.info(f"Processing sample: {sample_name}")
        threads = os.cpu_count() or 1
        bam = _open_alignment(alignment_file, reference_path=reference_path, threads=threads)
        logger.info(f"Opened alignment file with {threads} threads")
        extracted_data = []

        with Progress(transient=True) as progress:
            task = progress.add_task("Processing regions", total=len(repeats_df))
            for chrom, chrom_df in repeats_df.groupby("Chromosome", sort=False):
                logger.info(f"Processing chromosome: {chrom}")
                for _, row in chrom_df.iterrows():
                    start, end = int(row["Start"]), int(row["End"])
                    expected_repeat = row["Expected_Repeat"]
                    repeat_start, repeat_end = int(row["Repeat_Start"]), int(row["Repeat_End"])
                    logger.info(f"Processing region: {chrom}:{start}-{end}")
                    reads = bam.fetch(contig=chrom, start=start, stop=end)
                    for read in reads:
                        qualities = read.query_qualities if read.query_qualities else []
                        quality_stats = calculate_base_quality_stats(qualities)
                        read_start = read.reference_start
                        read_end = read.reference_start + read.query_length
                        umi = None
                        for tag in ("RX", "UR", "UB"):
                            try:
                                umi = read.get_tag(tag)
                                break
                            except KeyError:
                                continue
                        if umi is None:
                            umi = "NA"
                        if read_start <= repeat_start and read_end >= repeat_end:
                            if (
                                (
                                    quality_stats["Mean_Quality"] is not None
                                    and quality_stats["Mean_Quality"] > bq_threshold
                                )
                                and read.mapping_quality >= mq_threshold
                                and (keep_n or "N" not in read.query_sequence)
                            ):
                                extracted_data.append(
                                    {
                                        "Chromosome": chrom,
                                        "Region_Start": start,
                                        "Region_End": end,
                                        "Read_Start": read_start,
                                        "Read_End": read_end,
                                        "Read_Name": read.query_name,
                                        "Read_Sequence": read.query_sequence,
                                        "UMI": umi,
                                        "Mapping_Quality": read.mapping_quality,
                                        "Read_Length": (
                                            len(read.query_sequence) if read.query_sequence else 0
                                        ),
                                        "Expected_Repeat": expected_repeat,
                                        "Repeat_Coordinates": f"{chrom}:{repeat_start}-{repeat_end}",
                                        **quality_stats,
                                    }
                                )
                    progress.advance(task)

        extracted_df = pd.DataFrame(extracted_data)
        extracted_df["UMI"] = extracted_df["UMI"].fillna("NA")
        grouped = extracted_df.groupby(["Chromosome", "Repeat_Coordinates"])
        final_rows = []

        for (chrom, repeat_coord), group in track(
            grouped,
            total=grouped.ngroups,
            description="Deduplicating",
        ):
            # Collapse reads into UMI+position families and derive consensus allele
            dedup = collapse_umi_families(group)

            if len(dedup) >= min_reads:
                final_rows.append(dedup)
                continue

            all_reads = []
            reads = bam.fetch(
                contig=chrom,
                start=int(group["Region_Start"].iloc[0]),
                stop=int(group["Region_End"].iloc[0]),
            )
            for read in reads:
                read_start = read.reference_start
                read_end = read.reference_start + read.query_length
                if read_start > int(group["Read_Start"].iloc[0]) or read_end < int(
                    group["Read_End"].iloc[0]
                ):
                    continue
                qualities = read.query_qualities if read.query_qualities else []
                q_stats = calculate_base_quality_stats(qualities)
                if not keep_n and read.query_sequence and "N" in read.query_sequence:
                    continue
                umi = None
                for tag in ("RX", "UR", "UB"):
                    try:
                        umi = read.get_tag(tag)
                        break
                    except KeyError:
                        continue
                if umi is None:
                    umi = "NA"
                all_reads.append(
                    {
                        "Chromosome": chrom,
                        "Region_Start": int(group["Region_Start"].iloc[0]),
                        "Region_End": int(group["Region_End"].iloc[0]),
                        "Read_Start": read_start,
                        "Read_End": read_end,
                        "Read_Name": read.query_name,
                        "Read_Sequence": read.query_sequence,
                        "UMI": umi,
                        "Mapping_Quality": read.mapping_quality,
                        "Read_Length": len(read.query_sequence) if read.query_sequence else 0,
                        "Expected_Repeat": group["Expected_Repeat"].iloc[0],
                        "Repeat_Coordinates": repeat_coord,
                        **q_stats,
                    }
                )

            fallback_df = pd.DataFrame(all_reads)
            supplemental = pd.DataFrame()
            if not fallback_df.empty:
                fallback_df["UMI"] = fallback_df["UMI"].fillna("NA")
                # Collapse fallback reads into families as well
                fallback_df = collapse_umi_families(fallback_df)
                # Score and sort fallback consensus reads by quality
                fallback_df["Quality_Score"] = (
                    fallback_df["Mapping_Quality"] + fallback_df["Mean_Quality"]
                )
                fallback_df = fallback_df.sort_values(by="Quality_Score", ascending=False)

                # Calculate how many more reads needed
                n_missing = max(0, min_reads - len(dedup))

                # Take top-scoring consensus reads
                supplemental = fallback_df.head(n_missing)

            # Append deduplicated and supplemental reads for this locus
            combined = pd.concat([dedup, supplemental], ignore_index=True)
            final_rows.append(combined)

        final_df = pd.concat(final_rows, ignore_index=True)
        final_df.to_csv(output_file, index=False)
        bam.close()
        logger.info(f"Deduplicated reads saved to: {output_file}")

    except Exception as e:
        logger.error(f"An error occurred while processing the alignment file: {e}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract reads from an alignment file for specific repeat coordinates."
    )
    parser.add_argument(
        "-a", "--alignment", required=True, help="Path to the alignment file (.bam or .cram)."
    )
    parser.add_argument(
        "-r", "--repeats", required=True, help="Path to the CSV file with repeat coordinates."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to save the extracted reads CSV file."
    )
    parser.add_argument(
        "-g", "--reference", help="Path to the reference genome FASTA (required for CRAM input)."
    )
    parser.add_argument(
        "--bq_threshold", type=float, default=38, help="Base quality threshold (default: 38)"
    )
    parser.add_argument(
        "--mq_threshold", type=int, default=58, help="Mapping quality threshold (default: 58)"
    )
    parser.add_argument(
        "--keep_n",
        type=lambda x: (str(x).lower() == "false"),
        default=False,
        help="Whether to keep reads containing 'N' bases (default: False)",
    )
    parser.add_argument(
        "--min_reads",
        type=int,
        default=10,
        help="Minimum number of deduplicated reads per region (default: 10)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug-level logging.")
    parser.add_argument("--info", action="store_true", help="Enable info-level logging.")

    args = parser.parse_args()
    log_level = logging.DEBUG if args.verbose else logging.INFO if args.info else logging.WARNING
    logging.basicConfig(level=log_level, format="%(levelname)s:%(message)s")
    logger.setLevel(log_level)
    if args.alignment.lower().endswith(".cram") and not args.reference:
        raise ValueError("The --reference argument is required when the alignment input is CRAM.")
    repeat_coords = load_repeat_coordinates(args.repeats)
    repeat_coords = check_chr_format(args.alignment, repeat_coords, reference_path=args.reference)
    repeat_coords = repeat_coords.sort_values(by=["Chromosome", "Start"]).reset_index(drop=True)
    extract_reads_from_alignment(
        args.alignment,
        repeat_coords,
        args.output,
        args.bq_threshold,
        args.mq_threshold,
        args.keep_n,
        args.min_reads,
        reference_path=args.reference,
    )
