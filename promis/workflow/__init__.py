"""Helpers for locating PROMIS workflow assets."""

from __future__ import annotations

from importlib import resources
from pathlib import Path


def get_workflow_path() -> str:
    """Return the path to the installed workflow directory."""

    return str(resources.files(__name__))


def get_snakefile_path() -> str:
    """Return the absolute path to the packaged Snakefile."""

    return str(resources.files(__name__).joinpath("Snakefile"))


def get_default_config_path() -> str:
    """Return the absolute path to the default configuration file."""

    return str(resources.files(__name__).joinpath("config.yaml"))


def get_conda_env_path() -> str:
    """Return the absolute path to the rule-specific Conda environment file."""

    return str(resources.files(__name__).joinpath("environment.yml"))


def get_preset_path(name: str) -> str:
    """Return the packaged configuration template for ``name``."""

    return str(resources.files(__name__).joinpath("presets", f"{name}.yaml"))


def list_presets() -> list[str]:
    """Return the available packaged preset names."""

    return sorted(
        path.name.removesuffix(".yaml")
        for path in resources.files(__name__).joinpath("presets").iterdir()
        if path.name.endswith(".yaml")
    )


def resolve_resource_path(
    value: str | None,
    default_relative: str,
    packaged_default_path: str,
    run_dir: str | Path | None = None,
) -> str:
    """Resolve packaged defaults from the workflow and custom paths from the run directory."""

    selected = str(value or default_relative)
    candidate = Path(selected).expanduser()
    if candidate.is_absolute():
        return str(candidate)

    if candidate.as_posix() == Path(default_relative).as_posix():
        return str(Path(packaged_default_path).resolve())

    if candidate.parts and candidate.parts[0] == "database":
        workflow_root = Path(packaged_default_path).parents[1]
        packaged_candidate = workflow_root / candidate
        if packaged_candidate.exists():
            return str(packaged_candidate.resolve())

    base_dir = Path.cwd() if run_dir is None else Path(run_dir)
    return str((base_dir / candidate).resolve())


__all__ = [
    "get_workflow_path",
    "get_snakefile_path",
    "get_default_config_path",
    "get_conda_env_path",
    "get_preset_path",
    "list_presets",
    "resolve_resource_path",
]
