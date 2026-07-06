from __future__ import annotations

from pathlib import Path

import pytest

from promis import cli
from promis.workflow import get_snakefile_path


def test_copy_config_writes_default_config(tmp_path, capsys) -> None:
    config = tmp_path / "config.yaml"

    assert cli.main(["--copy-config", str(config)]) == 0

    text = config.read_text(encoding="utf-8")
    for key in ["output_dir", "alignment_files", "input_dir", "min_reads", "min_dev_percent"]:
        assert f"{key}:" in text
    assert "Wrote default PROMIS config" in capsys.readouterr().out


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


def test_cli_uses_launch_directory_as_default_workdir(tmp_path, monkeypatch) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("alignment_files: ''\noutput_dir: results/promis\n", encoding="utf-8")
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

    assert (
        cli.main(["--configfile", "config.yaml", "--dry-run", "--cores", "1", "--no-use-conda"])
        == 0
    )

    command = captured["command"]
    assert captured["cwd"] == str(tmp_path.resolve())
    assert command[command.index("--snakefile") + 1] == str(Path(get_snakefile_path()))
    assert command[command.index("--configfile") + 1] == str(config.resolve())
    assert "--use-conda" not in command
    assert captured["env"]["PROMIS_WORKFLOW_DIR"] == str(Path(get_snakefile_path()).parent)


def test_cli_respects_explicit_workdir(tmp_path, monkeypatch) -> None:
    project = tmp_path / "project"
    run_dir = tmp_path / "run"
    project.mkdir()
    run_dir.mkdir()
    config = project / "config.yaml"
    config.write_text("alignment_files: ''\noutput_dir: results/promis\n", encoding="utf-8")
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

    assert (
        cli.main(["--configfile", "config.yaml", "--workdir", str(run_dir), "--no-use-conda"]) == 0
    )

    assert captured["cwd"] == str(run_dir.resolve())
