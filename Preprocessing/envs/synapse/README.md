# Synapse Environment

Pattern: `requirements-complete`

Config variable: `SYNAPSE_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Preprocessing/envs/synapse/environment.yml \
  -p "${SYNAPSE_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${SYNAPSE_ENV}" python=3.10 pip
conda run -p "${SYNAPSE_ENV}" python -m pip install \
  -r Preprocessing/envs/synapse/requirements.txt
```

## Notes

- Authenticate after creation with `synapse login --auth-token <TOKEN>`.
