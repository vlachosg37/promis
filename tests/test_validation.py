from pathlib import Path

from promis.validation import build_sample_map, collect_alignment_files, validate_config


def test_build_sample_map_detects_duplicate_basenames() -> None:
    samples, duplicates = build_sample_map(["run1/sample.bam", "run2/sample.bam"])

    assert samples == {"sample": "run1/sample.bam"}
    assert duplicates == ["sample"]


def test_collect_alignment_files_from_input_dir_is_sorted(tmp_path) -> None:
    data = tmp_path / "data"
    data.mkdir()
    (data / "b.bam").write_text("bam\n", encoding="utf-8")
    (data / "a.cram").write_text("cram\n", encoding="utf-8")

    alignments = collect_alignment_files({"input_dir": "data"}, tmp_path)

    assert alignments == [
        str((data / "a.cram").resolve()),
        str((data / "b.bam").resolve()),
    ]


def test_collect_alignment_files_accepts_comma_and_space_separated_paths(tmp_path) -> None:
    first = tmp_path / "first.bam"
    second = tmp_path / "second.bam"

    alignments = collect_alignment_files({"alignment_files": f"{first}, {second}"}, tmp_path)

    assert alignments == [str(first), str(second)]


def test_validate_config_rejects_missing_index(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("not a real bam\n", encoding="utf-8")

    result = validate_config({"alignment_files": str(bam)}, run_dir=Path.cwd())

    assert not result.ok
    assert result.samples == {"sample": str(bam)}
    assert result.errors == [f"Alignment index not found: {bam}.bai"]


def test_validate_config_rejects_invalid_boolean(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("not a real bam\n", encoding="utf-8")

    result = validate_config({"alignment_files": str(bam), "use_GMM": "maybe"}, run_dir=Path.cwd())

    assert not result.ok
    assert any("use_GMM must be a boolean value" in error for error in result.errors)


def test_validate_config_rejects_conflicting_input_sources(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("bam\n", encoding="utf-8")
    (tmp_path / "sample.bam.bai").write_text("index\n", encoding="utf-8")

    result = validate_config({"alignment_files": [str(bam)], "input_dir": "data"}, run_dir=tmp_path)

    assert not result.ok
    assert "Set either alignment_files or input_dir, not both." in result.errors


def test_validate_config_rejects_ignored_legacy_keys(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("bam\n", encoding="utf-8")
    (tmp_path / "sample.bam.bai").write_text("index\n", encoding="utf-8")

    result = validate_config(
        {"alignment_files": [str(bam)], "collapse_umis": True}, run_dir=tmp_path
    )

    assert not result.ok
    assert "Unsupported configuration keys: collapse_umis" in result.errors
