# Stage 3 Integration Environment

Pattern: `hybrid`

Config variable: `BATCHCORR_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Processing/Tsai/Pipeline/envs/stage3_integration/environment.yml \
  -p "${BATCHCORR_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${BATCHCORR_ENV}" \
  -c conda-forge -c defaults \
  python=3.10 pip r-base=4.2 rpy2>=3.5
conda run -p "${BATCHCORR_ENV}" python -m pip install \
  -r Processing/Tsai/Pipeline/envs/stage3_integration/requirements.txt
```

## Notes

- ORA annotation loads the Mohammadi marker RDS via `rpy2`, so `r-base` must be available before installing the Python layer.
