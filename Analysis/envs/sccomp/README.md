# ACE sccomp Environment

Pattern: `hybrid`

Config variable: `SCCOMP_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Analysis/envs/sccomp/environment.yml \
  -p "${SCCOMP_ENV}"
conda run -p "${SCCOMP_ENV}" Rscript -e "install.packages('sccomp', repos='https://cloud.r-project.org')"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${SCCOMP_ENV}" \
  -c conda-forge -c defaults \
  r-base=4.2 r-dplyr r-tidyr r-readr r-ggplot2 r-patchwork pip
conda run -p "${SCCOMP_ENV}" Rscript -e "install.packages('sccomp', repos='https://cloud.r-project.org')"
conda run -p "${SCCOMP_ENV}" python -m pip install \
  -r Analysis/envs/sccomp/requirements.txt
```

## Notes

- `sccomp` is installed from CRAN because it is not bundled in the tracked conda spec.
