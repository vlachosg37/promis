from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

EXPECTED_RESULTS = {
    "toy_mss": {"score": 0.0, "unstable_regions": "0/5"},
    "toy_msi": {"score": 40.0, "unstable_regions": "2/5"},
}
COMBINED_COLUMNS = {
    "Score",
    "Unstable regions",
    "Total regions",
    "Call status",
    "Sample",
    "Evaluable_Loci",
    "Unstable_Loci",
    "Score_Percent",
    "Score_Fraction",
    "QC_Status",
}


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake is not installed")
def test_golden_workflow_scores_are_stable(tmp_path) -> None:
    pytest.importorskip("pysam")
    for sample in EXPECTED_RESULTS:
        if not Path(f"tests/data/{sample}.bam").exists():
            pytest.skip("golden toy BAM fixtures have not been generated")

    config_path = tmp_path / "config.golden.yaml"
    with Path("config/config.golden.yaml").open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    config["output_dir"] = str(tmp_path / "golden-promis")
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    subprocess.run(
        [
            "snakemake",
            "-s",
            "promis/workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "1",
        ],
        check=True,
    )

    combined = pd.read_csv(Path(config["output_dir"]) / "combined_results.csv")
    assert COMBINED_COLUMNS.issubset(combined.columns)
    for sample, expected in EXPECTED_RESULTS.items():
        row = combined.loc[combined["Sample"] == sample]
        assert len(row) == 1
        assert row.iloc[0]["Score"] == pytest.approx(expected["score"])
        assert row.iloc[0]["Unstable regions"] == expected["unstable_regions"]
        assert row.iloc[0]["Score_Percent"] == pytest.approx(expected["score"])
        assert row.iloc[0]["Score_Fraction"] == pytest.approx(expected["score"] / 100.0)
        assert row.iloc[0]["Evaluable_Loci"] == 5
        assert row.iloc[0]["QC_Status"] == "PASS"

        distribution = pd.read_csv(
            Path(config["output_dir"]) / sample / f"{sample}_distribution_analysis.csv"
        )
        callable_loci = distribution[distribution["Chromosome"] != "Summary"]
        summary_rows = distribution[distribution["Chromosome"] == "Summary"]
        assert len(callable_loci) == 5
        assert len(summary_rows) == 1

    resolved_config = Path(config["output_dir"]) / "resolved_config.yaml"
    run_metadata = Path(config["output_dir"]) / "run_metadata.json"
    assert resolved_config.exists()
    assert run_metadata.exists()
    metadata = json.loads(run_metadata.read_text(encoding="utf-8"))
    assert metadata["sample_count"] == 2
    assert metadata["thresholds"]["call_by"] == "both"
    assert set(metadata["samples"]) == set(EXPECTED_RESULTS)


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake is not installed")
def test_no_evaluable_loci_are_qc_failures(tmp_path) -> None:
    pytest.importorskip("pysam")
    if not Path("tests/data/toy_mss.bam").exists():
        pytest.skip("golden toy BAM fixtures have not been generated")

    config_path = tmp_path / "config.no-call.yaml"
    with Path("config/config.golden.yaml").open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    config["alignment_files"] = ["tests/data/toy_mss.bam"]
    config["output_dir"] = str(tmp_path / "no-call-promis")
    config["min_total_reads"] = 999
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    subprocess.run(
        [
            "snakemake",
            "-s",
            "promis/workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "1",
        ],
        check=True,
    )

    combined = pd.read_csv(Path(config["output_dir"]) / "combined_results.csv")
    row = combined.iloc[0]
    assert row["Sample"] == "toy_mss"
    assert row["Unstable regions"] == "0/0"
    assert row["Total regions"] == 0
    assert row["Evaluable_Loci"] == 0
    assert row["Unstable_Loci"] == 0
    assert pd.isna(row["Score"])
    assert pd.isna(row["Score_Percent"])
    assert pd.isna(row["Score_Fraction"])
    assert row["Call status"] == "NO_EVALUABLE_LOCI"
    assert row["QC_Status"] == "NO_EVALUABLE_LOCI"
