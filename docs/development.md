# Development and release

## Test from source

```bash
mamba env create -f envs/promis_test_env.yaml
mamba activate promis-test
python tests/create_tiny_bam.py
pytest
```

The committed tiny fixtures are synthetic. To verify the installed command-line
path, install the package into the environment and run a dry run from a new
directory:

```bash
python -m pip install --no-deps .
mkdir -p /tmp/promis_cli_test
cd /tmp/promis_cli_test
promis init --preset wes-wgs-hg38 --config config.yaml
```

## Release checklist

1. Merge the release changes and confirm CI is green on `main`.
2. Create the immutable GitHub release tag `v0.2.0`.
3. Download the tag archive and calculate its SHA-256.
4. Submit the matching `recipes/promis-msi/meta.yaml` to
   `bioconda/bioconda-recipes` and pass its lint, build, and isolated tests.
5. After Bioconda publishes the package, add the one-command Bioconda
   installation and package badge to the README.

The Bioconda recipe must test `promis --help`, `promis presets`, `promis init`,
`promis check`, `promis run --dry-run`, and `promis-find-ms-sites --help`.
