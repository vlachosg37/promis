from __future__ import annotations

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
    for sample, expected in EXPECTED_RESULTS.items():
        row = combined.loc[combined["Sample"] == sample]
        assert len(row) == 1
        assert row.iloc[0]["Score"] == pytest.approx(expected["score"])
        assert row.iloc[0]["Unstable regions"] == expected["unstable_regions"]

        distribution = pd.read_csv(
            Path(config["output_dir"]) / sample / f"{sample}_distribution_analysis.csv"
        )
        callable_loci = distribution[distribution["Chromosome"] != "Summary"]
        summary_rows = distribution[distribution["Chromosome"] == "Summary"]
        assert len(callable_loci) == 5
        assert len(summary_rows) == 1
