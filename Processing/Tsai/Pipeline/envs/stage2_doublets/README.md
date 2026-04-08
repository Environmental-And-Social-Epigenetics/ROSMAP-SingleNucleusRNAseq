# Stage 2 Doublet Environment

Pattern: `hybrid`

Config variable: `SINGLECELL_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Processing/Tsai/Pipeline/envs/stage2_doublets/environment.yml \
  -p "${SINGLECELL_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${SINGLECELL_ENV}" \
  -c conda-forge -c bioconda -c defaults \
  r-base=4.2 r-ggplot2 \
  bioconductor-biocparallel \
  bioconductor-genomeinfodbdata \
  bioconductor-singlecellexperiment \
  bioconductor-scdblfinder \
  bioconductor-zellkonverter
conda run -p "${SINGLECELL_ENV}" python -m pip install \
  -r Processing/Tsai/Pipeline/envs/stage2_doublets/requirements.txt
```
