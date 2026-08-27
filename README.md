# PROMIS: PROfiling of Microsatellite InStability

[![CI](https://github.com/vlachosg37/promis/actions/workflows/ci.yml/badge.svg)](https://github.com/vlachosg37/promis/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/vlachosg37/promis)](LICENSE)

PROMIS is a tumor-only, reference-free workflow for microsatellite instability
(MSI) analysis from WES, WGS, targeted-panel, and cfDNA sequencing data. It
reports a continuous MSI score together with evaluable-locus QC.

## Install

Bioconda installation will be added after `promis-msi` is accepted and
published. Until then, install from source:

```bash
git clone https://github.com/vlachosg37/promis.git
cd promis
mamba create -n promis -c conda-forge -c bioconda \
  python=3.12 pip pandas numpy scikit-learn pysam pyyaml rich tqdm numba snakemake
mamba activate promis
python -m pip install . --no-deps
```

## Run PROMIS

Create a project directory outside the source checkout, select a preset, edit
the generated configuration, then validate and run it:

```bash
mkdir promis_run
cd promis_run

promis presets
promis init --preset wes-wgs-hg38 --config config.yaml
nano config.yaml
promis check config.yaml
promis run config.yaml --cores 8
```

`alignment_files` accepts a YAML list or a comma-and-space-separated string.
Every BAM must be coordinate-sorted and indexed; CRAM input also requires the
matching reference FASTA and `.fai` index.

### Presets

| Preset | Use case |
| --- | --- |
| `wes-wgs-hg38` | Standard hg38 WES or WGS |
| `wes-wgs-cfdna-hg38` | hg38 cfDNA sequenced with WES/WGS-style loci and cfDNA thresholds |
| `qs-hrd-hg38` | QS-HRD targeted panel on hg38 |
| `tso500-hg19` | TSO500 targeted panel on hg19 |

The generated template documents all thresholds. For most analyses, edit only
the input paths, output directory, and any cohort-validated thresholds.

## Custom panels and locus discovery

Use `promis-find-ms-sites` to create a custom loci CSV. It scans all available
CPU cores by default; set `--cores` to the allocation on a shared node.

```bash
# Whole reference genome
promis-find-ms-sites --reference hg38.fa --output hg38_loci.csv --cores 8

# Restrict discovery to a panel BED file
promis-find-ms-sites --reference hg38.fa --output panel_loci.csv --bed panel.bed --cores 8

# Restrict discovery to regions with sufficient BAM coverage
promis-find-ms-sites \
  --reference hg38.fa --output covered_loci.csv --bam sample.bam \
  --min-coverage 30 --cores 8
```

Use the result in `config.yaml`:

```yaml
repeats: /absolute/path/to/custom_loci.csv
```

## Outputs

`results/promis/combined_results.csv` is the cohort summary. It begins with
`Sample`, followed by the MSI `Score`, region counts, call status, and QC
status. Per-sample read, repeat, and locus-level tables are written beneath the
same output directory.

## More documentation

- [Advanced run modes](docs/advanced-usage.md): direct Snakemake, Docker, and
  Apptainer/Singularity usage.
- [Development and release](docs/development.md): tests, release checks, and
  Bioconda submission workflow.

## Citation

Vlachos et al., *PROMIS: tumor-only profiling of microsatellite instability*,
bioRxiv (2025), DOI: TBD.
