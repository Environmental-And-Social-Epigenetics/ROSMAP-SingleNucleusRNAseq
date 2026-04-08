# Globus Environment

Pattern: `requirements-complete`

Config variable: `GLOBUS_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Preprocessing/envs/globus/environment.yml \
  -p "${GLOBUS_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${GLOBUS_ENV}" python=3.10 pip
conda run -p "${GLOBUS_ENV}" python -m pip install \
  -r Preprocessing/envs/globus/requirements.txt
```

## Notes

- Authenticate after creation with `globus login`.
