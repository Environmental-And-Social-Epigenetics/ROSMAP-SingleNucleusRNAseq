# COMPASS Environment

Pattern: `requirements-complete`

Config variable: `COMPASS_ANALYSIS_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Analysis/envs/compass/environment.yml \
  -p "${COMPASS_ANALYSIS_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${COMPASS_ANALYSIS_ENV}" python=3.10 pip
conda run -p "${COMPASS_ANALYSIS_ENV}" python -m pip install \
  -r Analysis/envs/compass/requirements.txt
```

## Notes

- Install IBM CPLEX separately and expose it to COMPASS before running analyses.
