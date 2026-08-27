import sys

import pytest

pysam = pytest.importorskip("pysam")

from promis.workflow.scripts.preprocess import find_MS_sites


def test_parse_args_accepts_explicit_core_count(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "promis-find-ms-sites",
            "--reference",
            "reference.fa",
            "--output",
            "loci.csv",
            "--cores",
            "3",
        ],
    )

    args = find_MS_sites.parse_args()

    assert args.cores == 3
