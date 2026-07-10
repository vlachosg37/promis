from __future__ import annotations

from pathlib import Path

import pytest

from promis import cli
from promis.workflow import get_snakefile_path


def _write_config(path: Path) -> None:
    path.write_text("alignment_files: ''\noutput_dir: results/promis\n", encoding="utf-8")


def _run_cli(tmp_path, monkeypatch, args):
    config = tmp_path / "config.yaml"
    _write_config(config)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        cli.shutil, "which", lambda name: "snakemake" if name == "snakemake" else None
    )

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


def test_copy_config_writes_default_config_in_empty_directory(
    tmp_path, monkeypatch, capsys
) -> None:
    monkeypatch.chdir(tmp_path)

    assert cli.main(["--copy-config", "config.yaml"]) == 0

    text = (tmp_path / "config.yaml").read_text(encoding="utf-8")
    for key in ["output_dir", "alignment_files", "input_dir", "min_reads", "min_dev_percent"]:
        assert f"{key}:" in text
    assert "Wrote default PROMIS config to config.yaml" in capsys.readouterr().out


def test_copy_config_refuses_overwrite_without_force(tmp_path) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("old: true\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc:
        cli.main(["--copy-config", str(config)])

    assert exc.value.code == 2
    assert config.read_text(encoding="utf-8") == "old: true\n"


def test_copy_config_overwrites_with_force(tmp_path) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("old: true\n", encoding="utf-8")

    assert cli.main(["--copy-config", str(config), "--force"]) == 0

    assert "output_dir:" in config.read_text(encoding="utf-8")


def test_print_config_still_prints_default_config(capsys) -> None:
    assert cli.main(["--print-config"]) == 0

    out = capsys.readouterr().out
    assert "output_dir:" in out
    assert "alignment_files:" in out


def test_cores_short_option_passes_snakemake_cores(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["-c", "8"])

    command = captured["command"]
    assert command[command.index("--cores") + 1] == "8"


def test_cores_long_option_accepts_all(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["--cores", "all"])

    command = captured["command"]
    assert command[command.index("--cores") + 1] == "all"


def test_jobs_short_option_is_passed_only_when_supplied(tmp_path, monkeypatch) -> None:
    without_jobs = _run_cli(tmp_path, monkeypatch, ["-c", "8"])
    assert "--jobs" not in without_jobs["command"]

    with_jobs = _run_cli(tmp_path, monkeypatch, ["-c", "8", "-j", "4"])
    command = with_jobs["command"]
    assert command[command.index("--jobs") + 1] == "4"


def test_jobs_long_option_is_passed(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["--jobs", "4"])

    command = captured["command"]
    assert command[command.index("--jobs") + 1] == "4"


def test_default_deployment_passes_no_software_deployment_method(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, [])

    command = captured["command"]
    assert "--use-conda" not in command
    assert "--software-deployment-method" not in command


@pytest.mark.parametrize(
    ("flag", "method"),
    [
        ("--use-conda", "conda"),
        ("--use-apptainer", "apptainer"),
        ("--use-singularity", "apptainer"),
    ],
)
def test_deployment_flags_pass_software_deployment_method(
    tmp_path, monkeypatch, flag, method
) -> None:
    captured = _run_cli(tmp_path, monkeypatch, [flag])

    command = captured["command"]
    assert command[command.index("--software-deployment-method") + 1] == method


def test_deployment_flags_are_mutually_exclusive(tmp_path) -> None:
    config = tmp_path / "config.yaml"
    _write_config(config)

    with pytest.raises(SystemExit) as exc:
        cli.main(["--configfile", str(config), "--use-conda", "--use-apptainer"])

    assert exc.value.code == 2


def test_cli_uses_launch_directory_as_default_workdir(tmp_path, monkeypatch) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["--dry-run", "--cores", "1"])

    command = captured["command"]
    config = tmp_path / "config.yaml"
    assert captured["cwd"] == str(tmp_path.resolve())
    assert command[command.index("--snakefile") + 1] == str(Path(get_snakefile_path()))
    assert command[command.index("--configfile") + 1] == str(config.resolve())
    assert captured["env"]["PROMIS_WORKFLOW_DIR"] == str(Path(get_snakefile_path()).parent)


def test_cli_respects_explicit_workdir(tmp_path, monkeypatch) -> None:
    project = tmp_path / "project"
    run_dir = tmp_path / "run"
    project.mkdir()
    run_dir.mkdir()
    config = project / "config.yaml"
    _write_config(config)
    monkeypatch.chdir(project)
    monkeypatch.setattr(
        cli.shutil, "which", lambda name: "snakemake" if name == "snakemake" else None
    )

    captured = {}

    def fake_run(command, cwd, env):
        captured["cwd"] = cwd

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(cli.subprocess, "run", fake_run)

    assert cli.main(["--configfile", "config.yaml", "--workdir", str(run_dir)]) == 0

    assert captured["cwd"] == str(run_dir.resolve())


def test_check_passes_for_existing_bam_and_index(tmp_path, capsys) -> None:
    bam = tmp_path / "sample.bam"
    bai = tmp_path / "sample.bam.bai"
    bam.write_text("not a real bam; validation only checks paths\n", encoding="utf-8")
    bai.write_text("index placeholder\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        f"alignment_files: {bam}\noutput_dir: results/promis\n",
        encoding="utf-8",
    )

    assert cli.main(["--configfile", str(config), "--check"]) == 0

    out = capsys.readouterr().out
    assert "PROMIS config check passed." in out
    assert "Samples: 1" in out
    assert "sample:" in out


def test_check_fails_for_duplicate_sample_names(tmp_path, capsys) -> None:
    run1 = tmp_path / "run1"
    run2 = tmp_path / "run2"
    run1.mkdir()
    run2.mkdir()
    for directory in (run1, run2):
        bam = directory / "sample.bam"
        bam.write_text("not a real bam\n", encoding="utf-8")
        (directory / "sample.bam.bai").write_text("index\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        f"alignment_files:\n  - {run1 / 'sample.bam'}\n  - {run2 / 'sample.bam'}\n",
        encoding="utf-8",
    )

    assert cli.main(["--configfile", str(config), "--check"]) == 1

    out = capsys.readouterr().out
    assert "ERROR: Duplicate sample names detected: sample" in out


def test_check_fails_for_cram_without_reference(tmp_path, capsys) -> None:
    cram = tmp_path / "sample.cram"
    cram.write_text("not a real cram\n", encoding="utf-8")
    (tmp_path / "sample.cram.crai").write_text("index\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(f"alignment_files: {cram}\n", encoding="utf-8")

    assert cli.main(["--configfile", str(config), "--check"]) == 1

    out = capsys.readouterr().out
    assert "ERROR: CRAM input requires reference_genome" in out


def test_init_writes_config_with_input_and_output_paths(tmp_path, capsys) -> None:
    config = tmp_path / "config.yaml"

    assert (
        cli.main(
            [
                "init",
                "--configfile",
                str(config),
                "--input-dir",
                "data",
                "--output-dir",
                "results/promis",
                "--mode",
                "wes",
            ]
        )
        == 0
    )

    text = config.read_text(encoding="utf-8")
    assert "input_dir: data" in text
    assert "output_dir: results/promis" in text
    assert "assay_mode: wes" in text
    assert "promis check --configfile" in capsys.readouterr().out


def test_cli_prints_run_summary(tmp_path, monkeypatch, capsys) -> None:
    captured = _run_cli(tmp_path, monkeypatch, ["--dry-run", "--cores", "2"])

    out = capsys.readouterr().out
    assert "PROMIS run" in out
    assert "Samples: 0" in out
    assert "Finished PROMIS run" in out
    assert captured["command"][captured["command"].index("--cores") + 1] == "2"
