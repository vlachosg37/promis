# PROMIS: PROfiling of Microsatellite InStability

[![CI](https://github.com/vlachosg37/promis/actions/workflows/ci.yml/badge.svg)](https://github.com/vlachosg37/promis/actions/workflows/ci.yml)

## Overview
PROMIS is a tumor-only, reference-free microsatellite instability (MSI) workflow built on Snakemake. It reports continuous MSI scores and evaluable-locus QC for whole-exome, whole-genome, targeted-panel, and cell-free DNA sequencing data.

## Features
- Tumor-only workflow (no matched normal required)
- Discrete mixture modeling at microsatellite loci
- Continuous MSI score with explicit evaluable-locus QC
- Reviewed presets for hg38 WES/WGS, QS-HRD, TSO500 hg19, and hg38 WES/WGS cfDNA
- Bundled loci and a reproducible conda environment

## Installation

### Bioconda

After the Bioconda recipe is merged, install PROMIS in one command:

```bash
mamba create -n promis -c conda-forge -c bioconda --strict-channel-priority promis-msi
```

### Local source install for development

```bash
git clone https://github.com/vlachosg37/promis.git
cd promis
mamba create -n promis-manual -c conda-forge -c bioconda \
  python=3.12 pip pandas numpy scikit-learn pysam \
  pyyaml rich tqdm numba snakemake-minimal
mamba activate promis-manual
python -m pip install . --no-deps
```

Run manual tests from a project directory outside the repository so relative
outputs and resources behave like an installed package.

## Running the pipeline

### Installed CLI mode

Run from a project directory, not from inside the installed package:

```bash
mkdir promis_run
cd promis_run

promis presets
promis init --preset wes-wgs-hg38 --config config.yaml
nano config.yaml
promis check config.yaml
promis run config.yaml --cores 8
```

By default, the CLI launches the packaged workflow and writes relative outputs
from the directory where you run `promis`. It does not enable conda or
container deployment unless you request it.

The `wes-wgs-cfdna-hg38` preset uses the hg38 WES/WGS loci with cfDNA-specific
thresholds for cfDNA produced with WES/WGS-style sequencing.

You can also create a pre-filled config for a project:

```bash
promis init --preset wes-wgs-hg38 --config config.yaml --input-dir /path/to/bams --output-dir results/promis
promis check config.yaml
promis run config.yaml --cores 8
```

Additional Snakemake options can be passed after `--`, for example `promis run config.yaml --cores 8 -- --profile slurm`.

For direct Snakemake usage with an installed package:

```bash
snakemake \
  -s "$(promis --workflow-dir)/Snakefile" \
  --configfile config.yaml \
  --cores 8
```

## Configuration
Place your `config.yaml` in the directory where you invoke `promis`, or pass it to `promis check` and `promis run`. Generate a complete template with `promis init --preset NAME --config config.yaml`; edit it to set:
- `output_dir`: destination for per-sample folders and combined results
- `alignment_files` **or** `input_dir`: explicit BAM/CRAM list or a directory to search recursively
- `repeats`, `scripts_dir`: override only if using custom resources
- Thresholds such as `min_reads`, `min_dev_reads`, `bq_threshold`, `mq_threshold`, `min_dev_percent`, and `use_GMM`

Comments in the template describe each option. Leave bundled resource paths such
as `repeats` and `scripts_dir` unchanged unless you have custom
resources; defaults resolve from the installed workflow, while custom relative
resource paths resolve from the directory where you run PROMIS. Most users only
need to edit `alignment_files`, `input_dir`, `output_dir`, `reference_genome`,
and thresholds.

Two portable config templates are provided:
- `config/config.example.yaml`: user-facing template with relative paths only.
- `config/config.test.yaml`: tiny synthetic fixture config for CI and local dry-runs.

## Inputs
- Coordinate-sorted, indexed BAM files for each tumor sample
- Reference genome FASTA with `.fai` only for CRAM input; BAM input does not need it
- MSI loci metadata provided in `promis/workflow/database/`

## Outputs
Results are organized under `output_dir` (default `results/promis`):
- `<sample>/<sample>_extracted_reads.csv`: filtered reads spanning each MSI locus
- `<sample>/<sample>_repeats_analysis.csv`: inferred repeat lengths per locus
- `<sample>/<sample>_distribution_analysis.csv`: stability calls and MSI status per locus
- `combined_results.csv`: cohort-level table summarizing MSI scores and unstable region counts
- `resolved_config.yaml` and `run_metadata.json`: fully resolved settings and provenance

## Testing

Create the test environment:

```bash
mamba env create -f envs/promis_test_env.yaml
mamba activate promis-test
```

Run local checks:

```bash
ruff check promis tests
black --check promis tests
python -m compileall promis tests
pytest
snakemake -s promis/workflow/Snakefile --configfile config/config.test.yaml -n --cores 1
snakemake -s promis/workflow/Snakefile --configfile config/config.golden.yaml --cores 1
```

The committed tiny fixture and generated golden fixtures are fully synthetic.
`tests/create_tiny_bam.py` creates 30 artificial reads over the 5 loci in
`tests/data/tiny_loci.csv` for each fixture sample; no patient, TCGA, WES, WGS,
panel, clinical, or local HPC data is used. The golden fixtures lock expected
algorithm behavior: `toy_mss` should score 0/5 unstable loci and `toy_msi`
should score 2/5 unstable loci.

The final PROMIS sample score is the percentage of unstable evaluable loci.
`combined_results.csv` preserves the legacy `Score` column and adds
`Score_Percent`, `Score_Fraction`, `Evaluable_Loci`, `Unstable_Loci`, and
`QC_Status` for unambiguous downstream use. `call_by: both` means the deviating
read count and deviating read percentage thresholds must both pass. The
`msi_deviation` setting is a read-level repeat-length shift threshold; it is not
the final sample score threshold. Samples with no evaluable loci are reported as
QC failures rather than true MSS calls.

Regenerate the synthetic BAM/BAI fixtures:

```bash
python tests/create_tiny_bam.py
```

Native Windows may not install `pysam` reliably. With Docker Desktop running:

```powershell
docker run --rm -v ${PWD}:/work -w /work python:3.12-slim sh -lc "pip install pysam && python tests/create_tiny_bam.py"
```

## Run modes and containers

PROMIS publishes Docker and native SIF images to GitHub Container Registry. Use
`latest` for quick testing; use version tags such as `v0.1.0` for real analyses
once the release tag exists.

The container includes only PROMIS code, runtime/test dependencies, and the tiny
fully synthetic test fixtures under `tests/data/`. No patient, TCGA, WES, WGS,
panel, clinical, or local HPC data is included.

### Cloned repository + Snakemake + conda

```bash
git clone https://github.com/vlachosg37/promis.git
cd promis

mamba env create -f promis/workflow/environment.yml
conda activate promis

cp promis/workflow/config.yaml config.yaml
nano config.yaml

snakemake \
  -s promis/workflow/Snakefile \
  --configfile config.yaml \
  --use-conda \
  --cores 8
```

### Cloned repository + Snakemake + Apptainer/Singularity

```bash
git clone https://github.com/vlachosg37/promis.git
cd promis

cp promis/workflow/config.yaml config.yaml
nano config.yaml

snakemake \
  -s promis/workflow/Snakefile \
  --configfile config.yaml \
  --software-deployment-method apptainer \
  --cores 8
```

Some clusters still call Apptainer `singularity`, and some Snakemake versions
may require `--software-deployment-method singularity`.

### Prebuilt SIF mode

```bash
singularity pull promis.sif oras://ghcr.io/vlachosg37/promis-sif:latest
singularity exec --bind "$PWD":"$PWD" --pwd "$PWD" promis.sif promis --copy-config config.yaml
nano config.yaml
```

Run with bind-mounted data and results:

```bash
singularity exec \
  --bind "$PWD":"$PWD" \
  --bind /path/to/data:/data \
  --pwd "$PWD" \
  promis.sif \
  promis --configfile config.yaml -c 8
```

Future pinned release form:

```bash
singularity pull promis_v0.1.0.sif oras://ghcr.io/vlachosg37/promis-sif:v0.1.0
```

### Docker mode

```bash
docker pull ghcr.io/vlachosg37/promis:latest
docker run --rm ghcr.io/vlachosg37/promis:latest promis --help
```

Run with bind-mounted project and data directories:

```bash
docker run --rm \
  -v /path/to/project:/project \
  -v /path/to/data:/data \
  -w /project \
  ghcr.io/vlachosg37/promis:latest \
  promis --configfile config.yaml -c 8
```

### Snakemake automatic Apptainer pull

Pin the image in your config when you need reproducible runs:

```yaml
container_image: "docker://ghcr.io/vlachosg37/promis:v0.1.0"
```

Then run through the installed CLI:

```bash
promis --configfile config.yaml -c 8 --use-apptainer
```

`--use-singularity` is accepted as an alias for clusters that still use that
name. Do not combine conda and Apptainer/Singularity deployment modes.

### CLI conda deployment

Use this mode when running outside a prebuilt container and you want Snakemake to
create rule-specific conda environments:

```bash
promis --configfile config.yaml -c 8 --use-conda
```

### Docker-to-Apptainer fallback

If the native SIF image is unavailable, pull the Docker image into a local SIF:

```bash
apptainer pull promis.sif docker://ghcr.io/vlachosg37/promis:latest
apptainer exec promis.sif \
  promis --configfile config.yaml -c 8
```

## Pre-release run-mode checklist

Installed CLI/source mode:

```bash
python -m pip install .
mkdir -p /tmp/promis_cli_test
cd /tmp/promis_cli_test
promis --copy-config config.yaml
promis --configfile config.yaml --dry-run -c 1
```

Repo Snakemake mode:

```bash
snakemake -s promis/workflow/Snakefile --configfile config/config.test.yaml -n --cores 1
snakemake -s promis/workflow/Snakefile --configfile config/config.test.yaml --cores 1
```

Container/SIF mode:

```bash
singularity pull --force promis.sif oras://ghcr.io/vlachosg37/promis-sif:latest
singularity exec promis.sif promis --help
singularity exec promis.sif promis --copy-config config.yaml --force
singularity exec promis.sif promis --configfile config.yaml --dry-run -c 1
```

Docker mode:

```bash
docker pull ghcr.io/vlachosg37/promis:latest
docker run --rm ghcr.io/vlachosg37/promis:latest promis --help
```

## Release

Do not create release tags until CI on `main` is green and the GHCR packages are
public and pullable. Version tags are preferred for reproducibility; `latest` is
only a convenience tag.

## Optional: Discovering microsatellite loci
PROMIS ships with a helper CLI, `promis-find-ms-sites`, to scan a reference genome for microsatellite loci:

```bash
promis-find-ms-sites \
  --reference hg38.fa \
  --output hg38_msi_loci.csv \
  --bam example.bam \
  --min-coverage 30
```

The resulting CSV can be used as a custom loci file in the PROMIS configuration.

## Citation
Vlachos et al., PROMIS: tumor-only profiling of microsatellite instability, bioRxiv (2025), DOI: TBD
