# ACE DEG / Nebula Environment

Pattern: `hybrid`

Config variable: `NEBULA_ENV`

## Conda setup

```bash
source config/paths.sh
conda env create \
  -f Analysis/envs/nebula/environment.yml \
  -p "${NEBULA_ENV}"
```

## Requirements setup

```bash
source config/paths.sh
conda create -y -p "${NEBULA_ENV}" \
  -c conda-forge -c bioconda -c defaults \
  python=3.10 pip r-base=4.2 r-dplyr bioconductor-zellkonverter \
  bioconductor-deseq2 bioconductor-singlecellexperiment bioconductor-scran \
  pandas>=2.0 numpy>=1.24 scipy>=1.10 h5py>=3.8
conda run -p "${NEBULA_ENV}" python -m pip install \
  -r Analysis/envs/nebula/requirements.txt
```

## Notes

- This environment supports the ACE DEG R workflow plus the Python H5AD split/prep helpers.
