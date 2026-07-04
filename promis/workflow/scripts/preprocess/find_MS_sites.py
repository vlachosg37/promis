"""
Microsatellite (Short Tandem Repeat) Loci Finder for Reference Genomes

This script scans a reference genome FASTA for all perfect microsatellite loci
(mono-, di-, tri-, and tetranucleotide repeats) across standard chromosomes
(chr1–22, chrX, chrY). It supports optional restriction of the search to
regions defined in a BED file or to regions of a BAM file that meet a minimum
coverage threshold. Results are written to a CSV compatible with the PROMIS
pipeline.

Key Features:
- Numba-accelerated search for perfect repeats of motifs 1–4 bp in length.
- User-defined minimum repeat thresholds for each motif size.
- Optional restriction of the search space using BED intervals or BAM coverage.
- Deduplication of overlapping loci, keeping the longest, leftmost instance per
  motif.
- Extraction of upstream and downstream sequence context for each locus.

Requirements:
- Python 3, pysam, pandas, numpy, numba, tqdm
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import pysam
from numba import njit
from tqdm import tqdm

# --- CONFIG ---
NUM_THREADS = 64
min_repeats = {1: 8, 2: 8, 3: 8, 4: 6}
max_motif = 4
context_length = 4  # bases for upstream/downstream context

standard_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

BASE_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Find microsatellites in a reference genome with optional restriction "
            "to BED regions or to regions of a BAM that meet a minimum coverage"
        )
    )
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA to scan")
    parser.add_argument("-o", "--output", required=True, help="Path to write the resulting CSV")
    parser.add_argument("--bed", help="Optional BED file restricting the search")
    parser.add_argument(
        "--bam",
        help="Optional BAM file. Only regions with at least --min-coverage are scanned",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=30,
        help="Minimum coverage required when --bam is provided",
    )
    args = parser.parse_args()
    if args.bed and args.bam:
        parser.error("Provide either --bed or --bam, not both")
    return args


def bed_intervals(bed_path):
    intervals = []
    with open(bed_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, *rest = line.rstrip().split("\t")
            intervals.append((chrom, int(start), int(end)))
    return intervals


def bam_coverage_intervals(bam_path, min_cov):
    intervals = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in bam.references:
            chrom_len = bam.get_reference_length(chrom)
            start = None
            for pileupcolumn in bam.pileup(chrom, 0, chrom_len, truncate=True):
                pos = pileupcolumn.pos
                cov = pileupcolumn.nsegments
                if cov >= min_cov:
                    if start is None:
                        start = pos
                else:
                    if start is not None:
                        intervals.append((chrom, start, pos))
                        start = None
            if start is not None:
                intervals.append((chrom, start, chrom_len))
    return intervals


def seq_to_array(seq):
    arr = np.empty(len(seq), dtype=np.uint8)
    for i, base in enumerate(seq):
        arr[i] = BASE_MAP.get(base, 4)  # non-ACGT as 4
    return arr


@njit
def find_str_loci(seq_arr, motif_lens, min_repeats_arr):
    results = []
    seq_len = seq_arr.size
    for motif_idx in range(len(motif_lens)):
        motif_len = motif_lens[motif_idx]
        repeat_min = min_repeats_arr[motif_idx]
        window = motif_len * repeat_min
        for i in range(seq_len - window + 1):
            motif = seq_arr[i : i + motif_len]
            # skip ambiguous motifs or mononucleotide polyN >1
            if np.any(motif > 3) or (motif_len > 1 and len(set(motif)) == 1):
                continue
            # count how many motif repeats
            n = 1
            while i + (n + 1) * motif_len <= seq_len:
                nextmotif = seq_arr[i + n * motif_len : i + (n + 1) * motif_len]
                if np.array_equal(nextmotif, motif):
                    n += 1
                else:
                    break
            if n >= repeat_min:
                repeat_len = motif_len * n
                results.append((i, i + repeat_len, motif_len, n))
    return results


def find_microsatellites(seq, chrom, interval_start, motif_sizes=(1, 2, 3, 4), min_repeats=None):
    motif_lens = np.array(list(motif_sizes), dtype=np.uint8)
    min_reps = np.array([min_repeats[m] for m in motif_lens], dtype=np.uint8)
    arr = seq_to_array(seq)
    loci = find_str_loci(arr, motif_lens, min_reps)
    result_list = []
    for start, end, motif_len, n in loci:
        motif = seq[start : start + motif_len]
        result_list.append(
            {
                "Chromosome": chrom,
                "Repeat_Start": interval_start + start,
                "Repeat_End": interval_start + end,
                "Motif": motif,
                "Motif_Length": motif_len,
                "Total_Length": end - start,
                "Num_Repeats": n,
                "Local_Start": start,
                "Local_End": end,
            }
        )
    return result_list


def deduplicate_repeats(repeats):
    df = pd.DataFrame(repeats)
    if df.empty:
        return df
    df["span"] = df["Repeat_End"] - df["Repeat_Start"]
    df = df.sort_values(
        ["Chromosome", "Motif", "Repeat_Start", "span"],
        ascending=[True, True, True, False],
    )
    result_rows = []
    seen_intervals = set()
    for _, row in df.iterrows():
        overlap = False
        for pos in range(row["Repeat_Start"], row["Repeat_End"]):
            if (row["Chromosome"], pos, row["Motif"]) in seen_intervals:
                overlap = True
                break
        if not overlap:
            result_rows.append(row)
            for pos in range(row["Repeat_Start"], row["Repeat_End"]):
                seen_intervals.add((row["Chromosome"], pos, row["Motif"]))
    return pd.DataFrame(result_rows)


def make_expected_repeat(motif, n):
    return f"({motif}){n}"


def expand_repeat(motif, n):
    return motif * n


def get_context_from_chromseq(chrom, start, end, ctx_len, chrom_seqs):
    seq = chrom_seqs[chrom]
    upstream = seq[max(0, start - ctx_len) : start]
    downstream = seq[end : end + ctx_len]
    return upstream, downstream


def process_interval(chrom, start, end, fasta_path, min_repeats, max_motif):
    fasta = pysam.FastaFile(fasta_path)
    seq = fasta.fetch(chrom, start, end)
    fasta.close()
    repeats = find_microsatellites(
        seq,
        chrom,
        start,
        motif_sizes=range(1, max_motif + 1),
        min_repeats=min_repeats,
    )
    for r in repeats:
        r["Sequence"] = seq[r["Local_Start"] : r["Local_End"]]
    return repeats


def main():
    args = parse_args()

    if args.bed:
        interval_list = bed_intervals(args.bed)
    elif args.bam:
        interval_list = bam_coverage_intervals(args.bam, args.min_coverage)
    else:
        interval_list = None

    with pysam.FastaFile(args.reference) as fasta:
        if interval_list is None:
            intervals_by_chrom = {
                chrom: [(0, fasta.get_reference_length(chrom))]
                for chrom in standard_chroms
                if chrom in fasta.references
            }
        else:
            intervals_by_chrom = {}
            for chrom, start, end in interval_list:
                if chrom in fasta.references:
                    intervals_by_chrom.setdefault(chrom, []).append((start, end))

    all_repeats = []
    tasks = []
    with ProcessPoolExecutor(max_workers=NUM_THREADS) as executor:
        for chrom, intervals in intervals_by_chrom.items():
            for start, end in intervals:
                tasks.append(
                    executor.submit(
                        process_interval, chrom, start, end, args.reference, min_repeats, max_motif
                    )
                )
        for fut in tqdm(as_completed(tasks), total=len(tasks), desc="Intervals", position=0):
            all_repeats.extend(fut.result())

    dedup_df = deduplicate_repeats(all_repeats)

    if not dedup_df.empty:
        dedup_df["Expected_Repeat"] = dedup_df.apply(
            lambda row: make_expected_repeat(row["Motif"], int(row["Num_Repeats"])),
            axis=1,
        )
        dedup_df["Expanded_Repeat"] = dedup_df.apply(
            lambda row: expand_repeat(row["Motif"], int(row["Num_Repeats"])),
            axis=1,
        )
        dedup_df["Start"] = dedup_df["Repeat_Start"]
        dedup_df["End"] = dedup_df["Repeat_End"]

        print("Loading full chromosome sequences into memory ...")
        chrom_seqs = {}
        with pysam.FastaFile(args.reference) as fasta:
            for chrom in tqdm(dedup_df["Chromosome"].unique(), desc="Chromosomes (seq load)"):
                chrom_seqs[chrom] = fasta.fetch(chrom)

        print("Starting context extraction using in-memory chromosome sequences ...")
        dedup_df[["Upstream_Context", "Downstream_Context"]] = dedup_df.apply(
            lambda row: pd.Series(
                get_context_from_chromseq(
                    row["Chromosome"],
                    int(row["Start"]),
                    int(row["End"]),
                    context_length,
                    chrom_seqs,
                )
            ),
            axis=1,
        )

        print("Context extraction done. Filtering loci without full context...")

        mask_up = dedup_df["Upstream_Context"].str.len() == context_length
        mask_down = dedup_df["Downstream_Context"].str.len() == context_length
        mask_both = mask_up & mask_down
        dedup_df = dedup_df[mask_both]

        final_cols = [
            "Chromosome",
            "Start",
            "End",
            "Sequence",
            "Expected_Repeat",
            "Repeat_Start",
            "Repeat_End",
            "Expanded_Repeat",
            "Upstream_Context",
            "Downstream_Context",
        ]
        dedup_df = dedup_df[final_cols]
        dedup_df.to_csv(args.output, index=False)
        print(
            f"Found {len(dedup_df)} deduplicated microsatellite loci. Results saved to: {args.output}"
        )
    else:
        print("No microsatellites found in the reference genome.")


if __name__ == "__main__":
    main()
