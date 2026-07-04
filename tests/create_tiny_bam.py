from __future__ import annotations

import csv
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
LOCI = DATA / "tiny_loci.csv"
BAM = DATA / "tiny.bam"
BAI = DATA / "tiny.bam.bai"
CHR1_LENGTH = 25_000_000
READS_PER_LOCUS = 6


def make_read(name: str, start: int, sequence: str) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.flag = 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = [(0, len(sequence))]
    read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    return read


def synthetic_sequence(row: dict[str, str], variant: int) -> str:
    sequence = row["Sequence"]
    if variant >= 4:
        repeat = row["Expanded_Repeat"]
        sequence = sequence.replace(repeat, repeat[:-3], 1)
    prefix = ("ACGT" * 20)[:50]
    suffix = ("TGCA" * 20)[:50]
    return prefix + sequence + suffix


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for path in (BAM, BAI):
        path.unlink(missing_ok=True)

    with LOCI.open(newline="") as fh:
        loci = list(csv.DictReader(fh))

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": CHR1_LENGTH}]}
    with pysam.AlignmentFile(BAM, "wb", header=header) as bam:
        for locus_index, row in enumerate(loci, start=1):
            region_start = int(row["Start"])
            for read_index in range(READS_PER_LOCUS):
                sequence = synthetic_sequence(row, read_index)
                reference_start = region_start - 51 + read_index
                read_name = f"tiny_locus{locus_index}_read{read_index + 1}"
                bam.write(make_read(read_name, reference_start, sequence))

    pysam.index(str(BAM))
    print(f"created {BAM} and {BAI}")


if __name__ == "__main__":
    main()
