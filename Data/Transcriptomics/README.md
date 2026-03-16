# Transcriptomics Data

Canonical local directory for all transcriptomics data. After cloning, these directories are empty (data files are too large for git). Populate them using the scripts in [`Data_Access/`](../../Data_Access/).

## Directory Structure

```
Transcriptomics/
├── Tsai/
│   ├── FASTQs/                # Raw FASTQ files (480 patients, ~9 TB)
│   ├── Cellranger_Output/     # Cell Ranger 8.0.0 count matrices
│   ├── Cellbender_Output/     # Ambient RNA-corrected matrices (478 samples)
│   └── Processing_Outputs/    # QC-filtered, doublet-removed, annotated h5ad files
├── DeJager/
│   ├── FASTQs/                # Raw FASTQ files (~127 libraries)
│   ├── Cellranger_Output/     # Cell Ranger count matrices (~47 samples)
│   ├── Cellbender_Output/     # Ambient RNA-corrected matrices (~131 samples)
│   └── Processing_Outputs/    # QC-filtered, doublet-removed, annotated h5ad files
└── copy_data.sbatch           # SLURM script to copy data from legacy scattered locations
```

## How to Populate

The fastest way to get started is to download the CellBender outputs and/or annotated objects from the Tsai Lab NAS:

```bash
source config/paths.sh

# Download CellBender outputs (input to Processing pipeline)
cd Data_Access/Transcriptomics/Tsai_Server/Tsai
DATA_TYPES=Cellbender_Output ./download_from_nas.sh

# Or download only the final annotated object (Stage 3 output)
./download_processing_outputs_from_nas.sh 03_Integrated
```

For raw FASTQs, see the [Data_Access README](../../Data_Access/README.md):
- **Tsai FASTQs**: Transfer from MIT Engaging via Globus
- **DeJager FASTQs**: Download from Synapse (`syn21438684`)

## Path Variables

All paths are defined in [`config/paths.sh`](../../config/paths.sh):

| Variable | Points to |
|----------|-----------|
| `$TSAI_FASTQS` | `Tsai/FASTQs/` |
| `$TSAI_CELLRANGER` | `Tsai/Cellranger_Output/` |
| `$TSAI_CELLBENDER` | `Tsai/Cellbender_Output/` |
| `$DEJAGER_FASTQS_DIR` | `DeJager/FASTQs/` |
| `$DEJAGER_CELLRANGER` | `DeJager/Cellranger_Output/` |
| `$DEJAGER_CELLBENDER` | `DeJager/Cellbender_Output/` |
