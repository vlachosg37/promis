"""Command-line interface for the PROMIS Snakemake workflow."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from . import __version__
from .validation import load_config, validate_config
from .workflow import (
    get_default_config_path,
    get_preset_path,
    get_snakefile_path,
    get_workflow_path,
    list_presets,
)

DEFAULT_CONFIG_FILENAME = "config.yaml"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PROMIS launches the packaged Snakemake MSI workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("command", nargs="?", help="init, check, run, or presets")
    parser.add_argument("config_path", nargs="?", help="Configuration YAML for check or run.")
    parser.add_argument("-c", "--cores", default=1, help="Cores available to Snakemake.")
    parser.add_argument("-j", "--jobs", default=None, help="Maximum concurrent Snakemake jobs.")
    parser.add_argument(
        "--configfile",
        "--config",
        dest="configfile",
        default=DEFAULT_CONFIG_FILENAME,
        help="Configuration YAML path.",
    )
    parser.add_argument("--preset", choices=list_presets(), default="wes-wgs-hg38")
    parser.add_argument("--workdir", default=None, help="Snakemake working directory.")
    deployment_group = parser.add_mutually_exclusive_group()
    deployment_group.add_argument("--use-conda", action="store_true")
    deployment_group.add_argument("--use-apptainer", action="store_true")
    deployment_group.add_argument("--use-singularity", action="store_true")
    parser.add_argument("--conda-prefix", default=None)
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("-p", "--printshellcmds", action="store_true")
    parser.add_argument("--print-config", action="store_true")
    parser.add_argument("--copy-config", metavar="PATH", default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--workflow-dir", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--input-dir", default=None)
    parser.add_argument("--alignment-files", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--mode", choices=["wes", "wgs", "panel", "cfdna"], default=None)
    parser.add_argument("-v", "--version", action="version", version=f"PROMIS {__version__}")
    return parser


def _destination(path: str) -> Path:
    destination = Path(path).expanduser()
    return destination if destination.is_absolute() else destination.resolve()


def _copy_template(source: Path, destination: Path, force: bool) -> None:
    if destination.exists() and not force:
        raise ValueError(f"Refusing to overwrite existing file: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")


def _print_check(configfile: Path, workdir: Path) -> int:
    result = validate_config(load_config(configfile), run_dir=workdir)
    print(f"PROMIS config check: {configfile}")
    print(f"Workdir: {workdir}")
    print(f"Samples: {len(result.samples)}")
    for sample, alignment in result.samples.items():
        print(f"  {sample}: {alignment}")
    for warning in result.warnings:
        print(f"WARNING: {warning}")
    for error in result.errors:
        print(f"ERROR: {error}")
    if result.ok:
        print("PROMIS config check passed.")
        return 0
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args, extra_args = parser.parse_known_args(argv)
    commands = {"init", "check", "run", "presets"}
    if args.command not in commands | {None}:
        extra_args = [args.command, *extra_args]
        args.command = None

    default_config = Path(get_default_config_path())
    if args.print_config:
        sys.stdout.write(default_config.read_text(encoding="utf-8"))
        return 0
    if args.workflow_dir:
        print(get_workflow_path())
        return 0
    if args.command == "presets":
        print("\n".join(list_presets()))
        return 0

    template = Path(get_preset_path(args.preset))
    if args.copy_config:
        try:
            _copy_template(template, _destination(args.copy_config), args.force)
        except ValueError as exc:
            parser.error(str(exc))
        print(f"Wrote PROMIS {args.preset} config to {args.copy_config}")
        return 0

    if args.command == "init":
        config_target = args.config_path or args.configfile
        destination = _destination(config_target)
        try:
            _copy_template(template, destination, args.force)
        except ValueError as exc:
            parser.error(str(exc))
        if (
            args.input_dir is not None
            or args.alignment_files is not None
            or args.output_dir is not None
        ):
            config = load_config(destination)
            if args.input_dir is not None:
                config["input_dir"] = args.input_dir
                config["alignment_files"] = []
            if args.alignment_files is not None:
                config["alignment_files"] = [
                    item.strip() for item in args.alignment_files.split(",") if item.strip()
                ]
                config["input_dir"] = ""
            if args.output_dir is not None:
                config["output_dir"] = args.output_dir
            destination.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        print(f"Wrote PROMIS {args.preset} config to {destination}")
        print(f"Check it with: promis check {destination}")
        print(f"Run it with: promis run {destination} --cores 8")
        return 0

    configfile = _destination(args.config_path or args.configfile)
    if not configfile.exists():
        parser.error(f"Configuration file not found: {configfile}")
    workdir = _destination(args.workdir) if args.workdir else Path.cwd().resolve()
    if not workdir.exists():
        parser.error(f"Working directory not found: {workdir}")

    if args.check or args.command == "check":
        return _print_check(configfile, workdir)

    check_result = validate_config(load_config(configfile), run_dir=workdir)
    if not check_result.ok:
        return _print_check(configfile, workdir)
    snakemake_executable = shutil.which("snakemake")
    if snakemake_executable is None:
        parser.error("The 'snakemake' executable was not found in the current environment.")

    command = [
        snakemake_executable,
        "--snakefile",
        get_snakefile_path(),
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
    command.extend(arg for arg in extra_args if arg != "--")

    config = load_config(configfile)
    print("PROMIS run")
    print(f"Preset: {config.get('preset', 'custom')}")
    print(f"Samples: {len(check_result.samples)}")
    print(f"Output: {config.get('output_dir', 'results/promis')}")
    print(f"Workdir: {workdir}")
    print(f"Cores: {args.cores}")
    env = os.environ.copy()
    env.setdefault("PROMIS_WORKFLOW_DIR", str(Path(get_snakefile_path()).parent))
    result = subprocess.run(command, cwd=str(workdir), env=env)
    if result.returncode == 0:
        output_dir = Path(config.get("output_dir", "results/promis"))
        print("Finished PROMIS run")
        print(f"Cohort summary: {output_dir / 'combined_results.csv'}")
        print(f"Per-sample tables: {output_dir}")
    return result.returncode


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
