"""Command-line interface for the PROMIS Snakemake workflow."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .validation import load_config, validate_config
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
        "command",
        nargs="?",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-c",
        "--cores",
        default=1,
        help="Total cores available to Snakemake. Use 'all' to use all available cores.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        default=None,
        help="Maximum concurrent Snakemake jobs, mainly useful with cluster/executor/profile modes.",
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
        help=(
            "Working directory from which Snakemake should execute the workflow. "
            "Defaults to the current directory; relative config and output paths "
            "resolve from this directory."
        ),
    )
    deployment_group = parser.add_mutually_exclusive_group()
    deployment_group.add_argument(
        "--use-conda",
        action="store_true",
        help="Use Snakemake conda deployment.",
    )
    deployment_group.add_argument(
        "--use-apptainer",
        action="store_true",
        help="Use Snakemake Apptainer deployment.",
    )
    deployment_group.add_argument(
        "--use-singularity",
        action="store_true",
        help="Alias for Snakemake Apptainer/Singularity deployment.",
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
        "--check",
        action="store_true",
        help="Validate the PROMIS config and input files, then exit.",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="Input directory to write into a config created by 'promis init'.",
    )
    parser.add_argument(
        "--alignment-files",
        type=str,
        default=None,
        help="Comma-separated BAM/CRAM paths to write into a config created by 'promis init'.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory to write into a config created by 'promis init'.",
    )
    parser.add_argument(
        "--mode",
        choices=["wes", "wgs", "panel", "cfdna"],
        default=None,
        help="Assay mode label recorded when creating a config with 'promis init'.",
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
    if args.command not in {None, "init", "check"}:
        extra_args = [args.command, *extra_args]
        args.command = None

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
            f"promis --configfile {display_destination} -c 8{os.linesep}"
        )
        return 0

    if args.command == "init":
        destination = Path(args.configfile).expanduser()
        if not destination.is_absolute():
            destination = destination.resolve()
        if destination.exists() and not args.force:
            parser.error(f"Refusing to overwrite existing file: {destination}")
        config = load_config(default_config)
        if args.input_dir is not None:
            config["input_dir"] = args.input_dir
            config["alignment_files"] = ""
        if args.alignment_files is not None:
            config["alignment_files"] = args.alignment_files
            config["input_dir"] = ""
        if args.output_dir is not None:
            config["output_dir"] = args.output_dir
        if args.mode is not None:
            config["assay_mode"] = args.mode
        destination.parent.mkdir(parents=True, exist_ok=True)
        import yaml

        destination.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        sys.stdout.write(
            f"Wrote PROMIS config to {destination}{os.linesep}"
            f"Check it with:{os.linesep}"
            f"promis check --configfile {destination}{os.linesep}"
            f"Run it with:{os.linesep}"
            f"promis --configfile {destination} -c 8{os.linesep}"
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

    if args.check or args.command == "check":
        try:
            config = load_config(configfile)
            result = validate_config(config, run_dir=workdir)
        except ValueError as exc:
            parser.error(str(exc))

        sys.stdout.write(f"PROMIS config check: {configfile}{os.linesep}")
        sys.stdout.write(f"Workdir: {workdir}{os.linesep}")
        sys.stdout.write(f"Samples: {len(result.samples)}{os.linesep}")
        for sample, alignment in result.samples.items():
            sys.stdout.write(f"  {sample}: {alignment}{os.linesep}")
        for warning in result.warnings:
            sys.stdout.write(f"WARNING: {warning}{os.linesep}")
        for error in result.errors:
            sys.stdout.write(f"ERROR: {error}{os.linesep}")
        if result.ok:
            sys.stdout.write("PROMIS config check passed." + os.linesep)
            return 0
        return 1

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

    if args.jobs:
        command.extend(["--jobs", str(args.jobs)])
    if args.use_conda:
        command.extend(["--software-deployment-method", "conda"])
    if args.use_apptainer or args.use_singularity:
        command.extend(["--software-deployment-method", "apptainer"])
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

    config = load_config(configfile)
    check_result = validate_config(config, run_dir=workdir)
    sys.stdout.write(
        f"PROMIS run{os.linesep}"
        f"Samples: {len(check_result.samples)}{os.linesep}"
        f"Output: {config.get('output_dir', 'results/promis')}{os.linesep}"
        f"Workdir: {workdir}{os.linesep}"
        f"Cores: {args.cores}{os.linesep}"
    )
    result = subprocess.run(command, cwd=str(workdir), env=env)
    if result.returncode == 0:
        output_dir = config.get("output_dir", "results/promis")
        sys.stdout.write(
            f"Finished PROMIS run{os.linesep}"
            f"Cohort summary: {Path(output_dir) / 'combined_results.csv'}{os.linesep}"
            f"Per-sample reports: {output_dir}{os.linesep}"
        )
    else:
        sys.stdout.write(
            "PROMIS run failed. Re-run with --printshellcmds for command details, "
            "or inspect the Snakemake error above." + os.linesep
        )
    return result.returncode


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    sys.exit(main())
