"""Validation helpers for PROMIS configuration and input files."""

from __future__ import annotations

from dataclasses import dataclass
from glob import glob
from pathlib import Path
from typing import Any

import yaml

from .workflow import get_workflow_path, resolve_resource_path

RESOURCE_DEFAULTS = {
    "repeats": "database/MSI_loci_hg38_coordinates_metadata_exonic_chr_rem_artefacts.csv",
    "cytoband": "database/cytoBand_hg38.txt",
    "scripts_dir": "scripts",
}
VALID_CALL_BY = {"count", "percent", "both"}
NUMERIC_THRESHOLDS = {
    "min_reads": int,
    "min_dev_reads": int,
    "bq_threshold": float,
    "mq_threshold": int,
    "msi_deviation": float,
    "min_dev_percent": float,
    "min_length_percent": float,
    "balance_tolerance": float,
    "min_total_reads": int,
}
BOOLEAN_KEYS = {"use_GMM", "filter_common_unstable"}
TRUE_VALUES = {"1", "true", "yes", "y"}
FALSE_VALUES = {"0", "false", "no", "n"}


@dataclass(frozen=True)
class ValidationResult:
    """Validation messages and discovered sample inputs."""

    samples: dict[str, str]
    errors: list[str]
    warnings: list[str]

    @property
    def ok(self) -> bool:
        return not self.errors


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as fh:
        loaded = yaml.safe_load(fh) or {}
    if not isinstance(loaded, dict):
        raise ValueError("Configuration file must contain a YAML mapping.")
    return loaded


def parse_bool(value: Any, key: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in {0, 1}:
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in TRUE_VALUES:
            return True
        if lowered in FALSE_VALUES:
            return False
    raise ValueError(f"{key} must be a boolean value: true/false, yes/no, or 1/0.")


def parse_int(value: Any, key: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be an integer.") from exc
    if parsed < 0:
        raise ValueError(f"{key} must be non-negative.")
    return parsed


def parse_float(value: Any, key: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be numeric.") from exc
    if parsed < 0:
        raise ValueError(f"{key} must be non-negative.")
    return parsed


def parse_call_by(value: Any) -> str:
    parsed = str(value).strip().lower()
    if parsed not in VALID_CALL_BY:
        raise ValueError("call_by must be one of: both, count, percent.")
    return parsed


def collect_alignment_files(config: dict[str, Any], run_dir: str | Path) -> list[str]:
    alignment_value = config.get("alignment_files", config.get("bam_files", ""))
    if isinstance(alignment_value, list):
        alignment_files = [str(entry) for entry in alignment_value if str(entry).strip()]
    elif isinstance(alignment_value, str) and alignment_value:
        alignment_files = [entry.strip() for entry in alignment_value.split(",") if entry.strip()]
    else:
        alignment_files = []

    if not alignment_files and config.get("input_dir"):
        search_dir = Path(run_dir) / str(config["input_dir"])
        bam_matches = glob(str(search_dir / "**" / "*.bam"), recursive=True)
        cram_matches = glob(str(search_dir / "**" / "*.cram"), recursive=True)
        alignment_files = sorted(bam_matches + cram_matches)

    return [
        str((Path(run_dir) / path).resolve()) if not Path(path).is_absolute() else path
        for path in alignment_files
    ]


def sample_name_from_alignment(path: str) -> str:
    return Path(path).name.removesuffix(".bam").removesuffix(".cram")


def build_sample_map(alignment_files: list[str]) -> tuple[dict[str, str], list[str]]:
    samples: dict[str, str] = {}
    duplicates: list[str] = []
    for alignment in alignment_files:
        sample = sample_name_from_alignment(alignment)
        if sample in samples:
            duplicates.append(sample)
            continue
        samples[sample] = alignment
    return samples, sorted(set(duplicates))


def _index_path_for_alignment(path: Path) -> Path:
    if path.suffix.lower() == ".bam":
        return path.with_suffix(path.suffix + ".bai")
    return path.with_suffix(path.suffix + ".crai")


def validate_config(config: dict[str, Any], run_dir: str | Path | None = None) -> ValidationResult:
    run_path = Path.cwd() if run_dir is None else Path(run_dir)
    errors: list[str] = []
    warnings: list[str] = []

    alignment_files = collect_alignment_files(config, run_path)
    samples, duplicates = build_sample_map(alignment_files)

    if not alignment_files:
        errors.append("No alignment files found. Set alignment_files or input_dir.")
    if duplicates:
        errors.append("Duplicate sample names detected: " + ", ".join(duplicates))

    reference_genome = str(config.get("reference_genome", "") or "")
    for alignment in alignment_files:
        alignment_path = Path(alignment)
        suffix = alignment_path.suffix.lower()
        if suffix not in {".bam", ".cram"}:
            errors.append(f"Unsupported alignment format: {alignment}")
            continue
        if not alignment_path.exists():
            errors.append(f"Alignment file not found: {alignment}")
            continue
        index_path = _index_path_for_alignment(alignment_path)
        if not index_path.exists():
            warnings.append(f"Alignment index not found: {index_path}")
        if suffix == ".cram" and not reference_genome:
            errors.append(f"CRAM input requires reference_genome: {alignment}")

    if reference_genome and not Path(reference_genome).exists():
        errors.append(f"Reference genome not found: {reference_genome}")

    workflow_dir = Path(get_workflow_path())
    for key, default_relative in RESOURCE_DEFAULTS.items():
        resource_path = Path(
            resolve_resource_path(
                config.get(key, default_relative),
                default_relative,
                str(workflow_dir / default_relative),
                run_dir=run_path,
            )
        )
        if not resource_path.exists():
            errors.append(f"Configured {key} path not found: {resource_path}")

    call_by = str(config.get("call_by", "both"))
    try:
        parse_call_by(call_by)
    except ValueError as exc:
        errors.append(str(exc))

    for key in BOOLEAN_KEYS:
        if key not in config:
            continue
        try:
            parse_bool(config[key], key)
        except ValueError as exc:
            errors.append(str(exc))

    for key, caster in NUMERIC_THRESHOLDS.items():
        if key not in config:
            continue
        try:
            if caster is int:
                parse_int(config[key], key)
            else:
                parse_float(config[key], key)
        except ValueError as exc:
            errors.append(str(exc))

    return ValidationResult(samples=samples, errors=errors, warnings=warnings)
