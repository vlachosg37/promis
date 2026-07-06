"""Command-line interface for the PROMIS Snakemake workflow."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .workflow import (
    get_default_config_path,
    get_snakefile_path,
    get_workflow_path,
)

DEFAULT_CONFIG_FILENAME = "config.yaml"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "PROMIS launches the packaged Snakemake workflow for microsatellite "
            "instability analysis."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-j",
        "--cores",
        type=int,
        default=1,
        help="Number of cores/jobs to use when running Snakemake.",
    )
    parser.add_argument(
        "--configfile",
        type=str,
        default=DEFAULT_CONFIG_FILENAME,
        help="Path to a Snakemake configuration YAML file.",
    )
    parser.add_argument(
        "--workdir",
        type=str,
        default=None,
        help="Working directory from which Snakemake should execute the workflow.",
    )
    parser.add_argument(
        "--use-conda",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable Snakemake's --use-conda flag to create rule-specific environments.",
    )
    parser.add_argument(
        "--conda-prefix",
        type=str,
        default=None,
        help="Optional directory used by Snakemake to store Conda environments.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Perform a dry run without executing rules (passes --dry-run to Snakemake).",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue independent jobs after failures (Snakemake's --keep-going).",
    )
    parser.add_argument(
        "-p",
        "--printshellcmds",
        action="store_true",
        help="Print shell commands that Snakemake executes.",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print the packaged default configuration file and exit.",
    )
    parser.add_argument(
        "--copy-config",
        type=str,
        default=None,
        metavar="PATH",
        help="Copy the packaged default configuration file to PATH and exit.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow --copy-config to overwrite an existing file.",
    )
    parser.add_argument(
        "--workflow-dir",
        action="store_true",
        help="Print the path to the installed workflow directory and exit.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"PROMIS {__version__}",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args, extra_args = parser.parse_known_args(argv)

    snakefile_path = Path(get_snakefile_path())
    default_config = Path(get_default_config_path())

    if args.print_config:
        sys.stdout.write(default_config.read_text())
        return 0

    if args.workflow_dir:
        sys.stdout.write(str(get_workflow_path()) + os.linesep)
        return 0

    if args.copy_config:
        display_destination = args.copy_config
        destination = Path(args.copy_config).expanduser()
        if not destination.is_absolute():
            destination = destination.resolve()
        if destination.exists() and not args.force:
            parser.error(f"Refusing to overwrite existing file: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(default_config.read_text(), encoding="utf-8")
        sys.stdout.write(
            f"Wrote default PROMIS config to {display_destination}{os.linesep}"
            f"Edit this file, then run:{os.linesep}"
            f"promis --configfile {display_destination} --cores 8{os.linesep}"
        )
        return 0

    if args.configfile:
        configfile = Path(args.configfile).expanduser()
        if not configfile.is_absolute():
            configfile = configfile.resolve()
    else:
        configfile = default_config

    if not configfile.exists():
        parser.error(f"Configuration file not found: {configfile}")

    if args.workdir is None:
        workdir = Path.cwd().resolve()
    else:
        workdir = Path(args.workdir).expanduser()
        if not workdir.is_absolute():
            workdir = workdir.resolve()
    if not workdir.exists():
        parser.error(f"Working directory not found: {workdir}")

    snakemake_executable = shutil.which("snakemake")
    if snakemake_executable is None:
        parser.error(
            "The 'snakemake' executable was not found. Install snakemake-minimal "
            "or snakemake in the current environment."
        )

    command = [
        snakemake_executable,
        "--snakefile",
        str(snakefile_path),
        "--cores",
        str(args.cores),
        "--configfile",
        str(configfile),
    ]

    if args.use_conda:
        command.append("--use-conda")
    if args.conda_prefix:
        command.extend(["--conda-prefix", args.conda_prefix])
    if args.dry_run:
        command.append("--dry-run")
    if args.keep_going:
        command.append("--keep-going")
    if args.printshellcmds:
        command.append("--printshellcmds")

    passthrough = [arg for arg in extra_args if arg != "--"]
    command.extend(passthrough)

    env = os.environ.copy()
    env.setdefault("PROMIS_WORKFLOW_DIR", str(snakefile_path.parent))

    result = subprocess.run(command, cwd=str(workdir), env=env)
    return result.returncode


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    sys.exit(main())
