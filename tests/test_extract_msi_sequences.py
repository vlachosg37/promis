import pandas as pd
import pytest

from promis.workflow.scripts.extract_MSI_sequences import (
    calculate_base_quality_stats,
    check_chr_format,
    collapse_umi_families,
)


def test_calculate_base_quality_stats() -> None:
    assert calculate_base_quality_stats([40, 42, 44]) == {
        "Mean_Quality": 42,
        "Median_Quality": 42,
    }
    assert calculate_base_quality_stats([]) == {"Mean_Quality": None, "Median_Quality": None}


def test_check_chr_format_keeps_matching_chr_prefix() -> None:
    repeats = pd.DataFrame({"Chromosome": ["chr1"]})

    converted = check_chr_format("tests/data/tiny.bam", repeats)

    assert converted["Chromosome"].tolist() == ["chr1"]


def test_check_chr_format_adds_chr_prefix_for_chr_bam() -> None:
    repeats = pd.DataFrame({"Chromosome": ["1"]})

    converted = check_chr_format("tests/data/tiny.bam", repeats)

    assert converted["Chromosome"].tolist() == ["chr1"]


def test_cram_requires_reference() -> None:
    with pytest.raises(ValueError, match="reference genome path is required"):
        check_chr_format("sample.cram", pd.DataFrame({"Chromosome": ["chr1"]}))


def test_collapse_umi_families_keeps_majority_best_read() -> None:
    reads = pd.DataFrame(
        [
            {
                "UMI": "A",
                "Read_Start": 1,
                "Read_End": 10,
                "Read_Sequence": "AAAA",
                "Mapping_Quality": 60,
                "Mean_Quality": 40,
                "Read_Length": 4,
                "Read_Name": "kept",
            },
            {
                "UMI": "A",
                "Read_Start": 1,
                "Read_End": 10,
                "Read_Sequence": "AAAA",
                "Mapping_Quality": 50,
                "Mean_Quality": 40,
                "Read_Length": 4,
                "Read_Name": "lower_mapq",
            },
            {
                "UMI": "A",
                "Read_Start": 1,
                "Read_End": 10,
                "Read_Sequence": "TTTT",
                "Mapping_Quality": 60,
                "Mean_Quality": 45,
                "Read_Length": 4,
                "Read_Name": "minority",
            },
        ]
    )

    collapsed = collapse_umi_families(reads)

    assert collapsed["Read_Name"].tolist() == ["kept"]
