from pathlib import Path

import pandas as pd
import pytest

pysam = pytest.importorskip("pysam")


def test_tiny_loci_fixture_has_five_rows() -> None:
    loci = Path("tests/data/tiny_loci.csv")

    assert loci.exists()
    assert len(pd.read_csv(loci)) == 5


def test_tiny_bam_fixture_is_valid_and_indexed() -> None:
    bam = Path("tests/data/tiny.bam")
    bai = Path("tests/data/tiny.bam.bai")

    assert bam.exists()
    assert bai.exists()
    assert pysam.quickcheck(str(bam)) in (None, "")

    with pysam.AlignmentFile(bam, "rb") as fh:
        assert fh.references == ("chr1",)
        assert fh.lengths[0] >= 25_000_000
        assert fh.has_index()
        reads = list(fh.fetch(until_eof=True))

    assert len(reads) == 30
    assert all(not read.is_secondary for read in reads)
    assert all(not read.is_supplementary for read in reads)
    assert all(not read.is_duplicate for read in reads)
    assert all(read.mapping_quality >= 60 for read in reads)


def test_golden_toy_bam_fixtures_are_valid_and_indexed() -> None:
    for sample in ("toy_mss", "toy_msi"):
        bam = Path(f"tests/data/{sample}.bam")
        bai = Path(f"tests/data/{sample}.bam.bai")
        if not bam.exists() or not bai.exists():
            pytest.skip("golden toy BAM fixtures have not been generated")

        assert bam.exists()
        assert bai.exists()
        assert pysam.quickcheck(str(bam)) in (None, "")

        with pysam.AlignmentFile(bam, "rb") as fh:
            assert fh.references == ("chr1",)
            assert fh.lengths[0] >= 25_000_000
            assert fh.has_index()
            reads = list(fh.fetch(until_eof=True))

        assert len(reads) == 30
        assert all(not read.is_secondary for read in reads)
        assert all(not read.is_supplementary for read in reads)
        assert all(not read.is_duplicate for read in reads)
        assert all(read.mapping_quality >= 60 for read in reads)
