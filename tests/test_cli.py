from __future__ import annotations

from pathlib import Path

import pytest

from promis import cli


def _write_config(path: Path, bam: Path) -> None:
    path.write_text(
        f"preset: wes-wgs-hg38\nalignment_files:\n  - {bam}\noutput_dir: results/promis\n",
        encoding="utf-8",
    )


def _run_cli(tmp_path, monkeypatch, args):
    bam = tmp_path / "sample.bam"
    bam.write_text("validation fixture\n", encoding="utf-8")
    (tmp_path / "sample.bam.bai").write_text("index\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    _write_config(config, bam)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cli.shutil, "which", lambda name: "snakemake" if name == "snakemake" else None)
    captured = {}

    def fake_run(command, cwd, env):
        captured["command"] = command
        captured["cwd"] = cwd
        captured["env"] = env

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    assert cli.main(["--configfile", "config.yaml", *args]) == 0
    return captured


def test_presets_lists_all_supported_templates(capsys) -> None:
    assert cli.main(["presets"]) == 0
    assert capsys.readouterr().out.splitlines() == [
        "qs-hrd-hg38",
        "tso500-hg19",
        "wes-wgs-cfdna-hg38",
        "wes-wgs-hg38",
    ]


def test_init_preserves_commented_cfdna_template(tmp_path, capsys) -> None:
    config = tmp_path / "cfdna.yaml"
    assert cli.main(["init", "--preset", "wes-wgs-cfdna-hg38", "--config", str(config)]) == 0
    text = config.read_text(encoding="utf-8")
    assert "cfDNA sequenced with WES/WGS-style data" in text
    assert "reference_genome: \"\"  # Required only for CRAM input" in text
    assert "min_dev_percent: 4.0" in text
    assert "promis run" in capsys.readouterr().out


def test_copy_config_refuses_overwrite_without_force(tmp_path) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("old: true\n", encoding="utf-8")
    with pytest.raises(SystemExit) as exc:
        cli.main(["--copy-config", str(config)])
    assert exc.value.code == 2


def test_check_passes_for_existing_bam_and_index(tmp_path, capsys) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("validation fixture\n", encoding="utf-8")
    (tmp_path / "sample.bam.bai").write_text("index\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    _write_config(config, bam)
    assert cli.main(["check", str(config)]) == 0
    assert "PROMIS config check passed." in capsys.readouterr().out


def test_check_rejects_no_input(tmp_path, capsys) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("preset: wes-wgs-hg38\nalignment_files: []\ninput_dir: ''\n", encoding="utf-8")
    assert cli.main(["check", str(config)]) == 1
    assert "No alignment files found" in capsys.readouterr().out


def test_run_subcommand_passes_cores(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["run", "config.yaml", "--cores", "8"])
    assert captured["command"][captured["command"].index("--cores") + 1] == "8"


def test_legacy_invocation_still_runs(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["-c", "2"])
    assert captured["command"][captured["command"].index("--cores") + 1] == "2"


def test_passthrough_options_follow_separator(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["run", "config.yaml", "--", "--profile", "slurm"])
    assert captured["command"][-2:] == ["--profile", "slurm"]
