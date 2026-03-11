# Data

External data files used by the pipeline and downstream analyses.

## Contents

| Directory | Description |
|-----------|-------------|
| `Phenotypes/` | Clinical phenotype data from ROSMAP (demographics, AD diagnosis, ACE scores, neuropathology, biomarkers) |

## Notes

- Raw sequencing data (FASTQs) is **not** stored here — it lives on cluster storage (see `config/paths.sh`)
- Preprocessed count matrices are stored at `${DEJAGER_PREPROCESSED}` and `${TSAI_PREPROCESSED}`
- Gene marker references are in `Processing/Tsai/Pipeline/Resources/`
