# Advanced run modes

The installed CLI in the project README is the recommended path. These modes
are available when an environment or cluster policy requires them.

## Direct Snakemake

```bash
snakemake \
  -s "$(promis --workflow-dir)/Snakefile" \
  --configfile config.yaml \
  --cores 8
```

For rule-specific Conda environments, add `--use-conda` to `promis run`:

```bash
promis run config.yaml --cores 8 --use-conda
```

## Apptainer or Singularity

```bash
promis run config.yaml --cores 8 --use-apptainer
```

`--use-singularity` is an alias for clusters that use the older command name.
Do not combine Conda and Apptainer/Singularity deployment modes.

To run the published container directly, bind both the project and data
directories:

```bash
apptainer exec \
  --bind "$PWD":"$PWD" \
  --bind /path/to/data:/data \
  --pwd "$PWD" \
  docker://ghcr.io/vlachosg37/promis:latest \
  promis run config.yaml --cores 8
```

For reproducible analyses, use a versioned image tag rather than `latest`.

## Docker

```bash
docker run --rm \
  -v /path/to/project:/project \
  -v /path/to/data:/data \
  -w /project \
  ghcr.io/vlachosg37/promis:latest \
  promis run config.yaml --cores 8
```
