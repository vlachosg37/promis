# Test Data

`tests/create_tiny_bam.py` generates fully synthetic BAM/BAI fixtures with
30 artificial reads over 5 loci from `tests/data/tiny_loci.csv`.

- `tiny.bam`: smoke-test fixture with instability-like reads at each locus.
- `toy_mss.bam`: golden stable fixture expected to score 0/5 unstable loci.
- `toy_msi.bam`: golden mixed fixture expected to score 2/5 unstable loci.

These files contain no patient, TCGA, WES, WGS, panel, clinical, or local HPC
data.

To regenerate them:

```bash
conda env create -f envs/promis_test_env.yaml
conda activate promis-test
python tests/create_tiny_bam.py
```

Native Windows may not support `pysam` reliably. Use WSL, Ubuntu, or CI to
regenerate the BAM/BAI if conda cannot install `pysam` locally.

With Docker Desktop running on Windows:

```powershell
docker run --rm -v ${PWD}:/work -w /work python:3.12-slim sh -lc "pip install pysam && python tests/create_tiny_bam.py"
```
