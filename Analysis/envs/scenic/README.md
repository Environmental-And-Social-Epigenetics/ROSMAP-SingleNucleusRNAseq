# SCENIC Environment

Pattern: `requirements-complete`

Config variable: `SCENIC_ANALYSIS_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Analysis/envs/scenic/environment.yml \
  -p "${SCENIC_ANALYSIS_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${SCENIC_ANALYSIS_ENV}" python=3.10 pip
conda run -p "${SCENIC_ANALYSIS_ENV}" python -m pip install \
  -r Analysis/envs/scenic/requirements.txt
```

## Notes

- Download the motif ranking reference files separately before running pySCENIC.
