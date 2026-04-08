# bcftools Environment

Pattern: `hybrid`

Config variable: `BCFTOOLS_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Preprocessing/envs/bcftools/environment.yml \
  -p "${BCFTOOLS_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${BCFTOOLS_ENV}" \
  -c bioconda -c conda-forge -c defaults \
  bcftools>=1.17 samtools>=1.17 htslib>=1.17
conda run -p "${BCFTOOLS_ENV}" python -m pip install \
  -r Preprocessing/envs/bcftools/requirements.txt
```

## Notes

- This environment is intentionally conda-first because the core tools are not pip packages.
