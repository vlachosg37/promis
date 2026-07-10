"""PROMIS Snakemake pipeline packaging utilities."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("promis")
except PackageNotFoundError:  # pragma: no cover - fallback for local execution
    __version__ = "0.1.1"

from .workflow import (  # noqa: E402
    get_default_config_path,
    get_snakefile_path,
    get_workflow_path,
)

__all__ = [
    "__version__",
    "get_default_config_path",
    "get_snakefile_path",
    "get_workflow_path",
]
