"""Helpers for locating PROMIS workflow assets."""

from __future__ import annotations

from importlib import resources


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


__all__ = [
    "get_workflow_path",
    "get_snakefile_path",
    "get_default_config_path",
    "get_conda_env_path",
]
