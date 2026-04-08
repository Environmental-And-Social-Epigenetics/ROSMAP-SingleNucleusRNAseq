# Legacy Data Preparation Scripts

**Status: DEPRECATED**

These scripts were the original ad-hoc preprocessing workflow used before the
standardized Processing pipeline (`Processing/*/Pipeline/`) existed. They are
preserved for historical reference only.

**You do NOT need these scripts** if you have the annotated h5ad from the
Processing pipeline (`tsai_annotated.h5ad` / `dejager_annotated.h5ad`). The
current Processing pipeline handles all QC filtering, doublet removal,
integration, and annotation.

Scripts in this directory contain hardcoded paths to the decommissioned MIT
Openmind cluster and are not portable without modification.
