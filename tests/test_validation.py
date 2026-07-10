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


def test_validate_config_warns_for_missing_index(tmp_path) -> None:
    bam = tmp_path / "sample.bam"
    bam.write_text("not a real bam\n", encoding="utf-8")

    result = validate_config({"alignment_files": str(bam)}, run_dir=Path.cwd())

    assert result.ok
    assert result.samples == {"sample": str(bam)}
    assert result.warnings == [f"Alignment index not found: {bam}.bai"]
