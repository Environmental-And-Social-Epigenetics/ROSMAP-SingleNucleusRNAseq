# DEG Analysis Environment

Pattern: `hybrid`

Config variable: `DEG_ANALYSIS_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Analysis/envs/deg/environment.yml \
  -p "${DEG_ANALYSIS_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${DEG_ANALYSIS_ENV}" \
  -c conda-forge -c bioconda -c defaults \
  python=3.10 pip r-base=4.2 \
  bioconductor-deseq2 bioconductor-edger bioconductor-limma \
  bioconductor-biocparallel r-ggplot2 r-dplyr r-tidyr r-pheatmap
conda run -p "${DEG_ANALYSIS_ENV}" python -m pip install \
  -r Analysis/envs/deg/requirements.txt
```
