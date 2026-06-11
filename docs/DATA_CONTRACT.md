# Data Contract

`Data/Transcriptomics/{Tsai,DeJager}` is the canonical data namespace for this repository.

Two local modes are supported:

1. Copy data into the namespace.
2. Symlink directories inside the namespace to external storage, including `/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data/Single_Nucleus/`.

Tracked files in `Data/` are limited to READMEs, manifests, small phenotype/reference metadata, and small fixtures. FASTQ, BAM, VCF, H5, H5AD, RDS, MTX, Zarr, logs, reports, and run outputs are not tracked unless explicitly whitelisted.

`Data/manifest.tsv` records expected roots, patterns, approximate counts, required status, provenance, and notes.

