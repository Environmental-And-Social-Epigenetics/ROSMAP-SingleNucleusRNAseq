# Environment Specifications

This directory is the canonical environment inventory for the repository.

- `envs/preprocessing/`: data download, Cell Ranger helpers, CellBender, demuxlet support.
- `envs/processing/`: Stage 1 QC, Stage 2 doublet removal, Stage 3 integration/annotation.
- `envs/analysis/`: ACE and downstream analysis environments.

Create environments with `setup/install_envs.sh`. Historical env directories under `Preprocessing/`, `Processing/`, and `Analysis/` are compatibility references only.

