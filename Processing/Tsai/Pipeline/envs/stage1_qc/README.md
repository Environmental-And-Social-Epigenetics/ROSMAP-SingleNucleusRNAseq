# Stage 1 QC Environment

Pattern: `requirements-complete`

Config variable: `QC_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Processing/Tsai/Pipeline/envs/stage1_qc/environment.yml \
  -p "${QC_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${QC_ENV}" python=3.10 pip
conda run -p "${QC_ENV}" python -m pip install \
  -r Processing/Tsai/Pipeline/envs/stage1_qc/requirements.txt
```
