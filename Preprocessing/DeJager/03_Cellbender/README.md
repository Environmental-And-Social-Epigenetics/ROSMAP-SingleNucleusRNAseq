# 03_Cellbender

Scripted CellBender submission for the DeJager cohort.

## Entry Points

- `Run_DeJager_Cellbender.py`: discovers libraries and submits jobs
- `Run_DeJager_Cellbender.sh`: wrapper that sources `config/paths.sh`
- `example_cellbender.sh`: reference single-library script

## Quick Start

```bash
source config/paths.sh
bash setup/install_envs.sh --preprocessing

cd Preprocessing/DeJager/03_Cellbender
bash Run_DeJager_Cellbender.sh --dry-run
bash Run_DeJager_Cellbender.sh --submit
```

Useful options:

- `--library-ids LIB1 LIB2`
- `--limit N`
- `--overwrite`

## Inputs And Outputs

Input:

```text
${DEJAGER_COUNTS}/{LibraryID}/outs/raw_feature_bc_matrix.h5
```

Output:

```text
${DEJAGER_PREPROCESSED}/{LibraryID}/
  processed_feature_bc_matrix.h5
  processed_feature_bc_matrix_filtered.h5
  processed_feature_bc_matrix.pdf
```

Submission logs are written alongside the preprocessing output tree rather than
into the repo.
