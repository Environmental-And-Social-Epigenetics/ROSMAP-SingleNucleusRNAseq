# Data

Canonical local location for all pipeline data. After cloning the repository, only `Phenotypes/` is available immediately (small CSVs tracked in git). All transcriptomics data must be populated using the scripts in [`Data_Access/`](../Data_Access/).

## Contents

| Directory | Description | Git-tracked? |
|-----------|-------------|:------------:|
| `Phenotypes/` | Clinical phenotype data from ROSMAP (demographics, AD diagnosis, ACE scores, neuropathology, biomarkers) | Yes |
| `Transcriptomics/Tsai/FASTQs/` | Raw FASTQ files for 480 Tsai patients (~9 TB) | No |
| `Transcriptomics/Tsai/Cellranger_Output/` | Cell Ranger 8.0.0 count matrices (480 samples) | No |
| `Transcriptomics/Tsai/Cellbender_Output/` | Ambient RNA-corrected matrices (478 samples) | No |
| `Transcriptomics/Tsai/Processing_Outputs/` | QC-filtered, doublet-removed, and annotated h5ad files | No |
| `Transcriptomics/DeJager/FASTQs/` | Raw FASTQ files for ~127 DeJager libraries | No |
| `Transcriptomics/DeJager/Cellranger_Output/` | Cell Ranger count matrices (~47 samples) | No |
| `Transcriptomics/DeJager/Cellbender_Output/` | Ambient RNA-corrected matrices (~131 samples) | No |
| `Transcriptomics/DeJager/Processing_Outputs/` | QC-filtered, doublet-removed, and annotated h5ad files | No |

## Populating Data After Cloning

Transcriptomics data is too large for git. Use the scripts in `Data_Access/` to download from the Tsai Lab NAS, transfer via Globus, or download from Synapse:

```bash
source config/paths.sh

# Example: download CellBender outputs from NAS
cd Data_Access/Transcriptomics/Tsai_Server/Tsai
DATA_TYPES=Cellbender_Output ./download_from_nas.sh

# Example: download only the annotated h5ad
./download_processing_outputs_from_nas.sh 03_Integrated
```

See [`Data_Access/README.md`](../Data_Access/README.md) for the full reference of available data, storage locations, and transfer scripts.

## Notes

- All path variables (e.g., `$TSAI_FASTQS`, `$TSAI_CELLBENDER`) are defined in [`config/paths.sh`](../config/paths.sh)
- Gene marker references are in `Processing/Tsai/Pipeline/Resources/`
