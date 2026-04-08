# GSEA Environment

Pattern: `hybrid`

Config variable: `GSEA_ANALYSIS_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Analysis/envs/gsea/environment.yml \
  -p "${GSEA_ANALYSIS_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${GSEA_ANALYSIS_ENV}" \
  -c conda-forge -c bioconda -c defaults \
  r-base=4.2 r-ggplot2 r-dplyr \
  bioconductor-clusterprofiler \
  bioconductor-org.hs.eg.db \
  r-webgestaltr
conda run -p "${GSEA_ANALYSIS_ENV}" python -m pip install \
  -r Analysis/envs/gsea/requirements.txt
```
