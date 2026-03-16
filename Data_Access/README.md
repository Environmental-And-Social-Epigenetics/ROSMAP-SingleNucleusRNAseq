# Data Access

This directory contains scripts for transferring ROSMAP snRNA-seq data between storage locations, and serves as the central reference for where each dataset lives and how to access it.

All path variables referenced below (e.g., `$TSAI_FASTQS`) are defined in [`config/paths.sh`](../config/paths.sh). Source it before using any scripts:

```bash
source config/paths.sh
```

---

## Quick Reference: Where Is My Data?

| Data | Primary Location | Backup | How to Get It |
|------|-----------------|--------|---------------|
| **Tsai FASTQs** | Engaging cluster | NAS, Openmind | Globus: [`send_globus.sh`](Transcriptomics/Engaging-Openmind_Transfer/Tsai/send_globus.sh) |
| **DeJager FASTQs** | Synapse (`syn21438684`) | NAS, Openmind | Synapse: [`Preprocessing/DeJager/01_FASTQ_Download/`](../Preprocessing/DeJager/01_FASTQ_Download/) |
| **Tsai CellRanger** | Openmind (`$TSAI_CELLRANGER`) | NAS | NAS: [`download_from_nas.sh`](Transcriptomics/Tsai_Server/Tsai/download_from_nas.sh) |
| **DeJager CellRanger** | Openmind (`$DEJAGER_CELLRANGER`) | NAS | NAS: [`download_from_nas.sh`](Transcriptomics/Tsai_Server/DeJager/download_from_nas.sh) |
| **Tsai CellBender** | Openmind (`$TSAI_CELLBENDER`) | NAS | NAS: [`download_from_nas.sh`](Transcriptomics/Tsai_Server/Tsai/download_from_nas.sh) |
| **DeJager CellBender** | Openmind permanent (`$DEJAGER_CELLBENDER`) | NAS | NAS: [`download_from_nas.sh`](Transcriptomics/Tsai_Server/DeJager/download_from_nas.sh) |
| **Tsai QC-Filtered** | `$TSAI_QC_FILTERED` | NAS | NAS: [`download_processing_outputs_from_nas.sh`](Transcriptomics/Tsai_Server/Tsai/download_processing_outputs_from_nas.sh) or run pipeline Stage 1 |
| **Tsai Doublet-Removed** | `$TSAI_DOUBLET_REMOVED` | NAS | NAS: [`download_processing_outputs_from_nas.sh`](Transcriptomics/Tsai_Server/Tsai/download_processing_outputs_from_nas.sh) or run pipeline Stage 2 |
| **Tsai Annotated** | `$TSAI_INTEGRATED` | NAS | NAS: [`download_processing_outputs_from_nas.sh 03_Integrated`](Transcriptomics/Tsai_Server/Tsai/download_processing_outputs_from_nas.sh) or run pipeline Stage 3 |
| **DeJager QC-Filtered** | `$DEJAGER_QC_FILTERED` | NAS | NAS: [`download_processing_outputs_from_nas.sh`](Transcriptomics/Tsai_Server/DeJager/download_processing_outputs_from_nas.sh) or run pipeline Stage 1 |
| **DeJager Doublet-Removed** | `$DEJAGER_DOUBLET_REMOVED` | NAS | NAS: [`download_processing_outputs_from_nas.sh`](Transcriptomics/Tsai_Server/DeJager/download_processing_outputs_from_nas.sh) or run pipeline Stage 2 |
| **DeJager Annotated** | `$DEJAGER_INTEGRATED` | NAS | NAS: [`download_processing_outputs_from_nas.sh 03_Integrated`](Transcriptomics/Tsai_Server/DeJager/download_processing_outputs_from_nas.sh) or run pipeline Stage 3 |
| **Phenotypes** | `Data/Phenotypes/` (in repo) | NAS | Already in repo |

---

## Transcriptomics -- Raw Data

### Tsai FASTQs

| | |
|---|---|
| **Description** | 480 patients, 5,197 FASTQ files, ~9 TB |
| **Primary location** | Engaging cluster (scattered across `/nfs/picower*` filesystems) |
| **FASTQ index CSV** | `$TSAI_FASTQS_CSV` -- maps each FASTQ to its patient, library, and path on Engaging |
| **Canonical path (Openmind)** | `$TSAI_FASTQS` = `Data/Transcriptomics/Tsai/FASTQs/` |
| **NAS backup** | `LabShared/mabdel03/ROSMAP/Data/Transcriptomics/Tsai/FASTQs/` |

**How to get Tsai FASTQs:**

| Scenario | Script |
|----------|--------|
| Engaging -> Openmind | [`Engaging-Openmind_Transfer/Tsai/send_globus.sh`](Transcriptomics/Engaging-Openmind_Transfer/Tsai/send_globus.sh) |
| Openmind -> Engaging | [`Engaging-Openmind_Transfer/Tsai/receive_globus.sh`](Transcriptomics/Engaging-Openmind_Transfer/Tsai/receive_globus.sh) |
| NAS -> local | [`Tsai_Server/Tsai/download_from_nas.sh`](Transcriptomics/Tsai_Server/Tsai/download_from_nas.sh) |
| Local -> NAS | [`Tsai_Server/Tsai/upload_to_nas.sh`](Transcriptomics/Tsai_Server/Tsai/upload_to_nas.sh) (parallel: [`slurm_upload.sh`](Transcriptomics/Tsai_Server/Tsai/slurm_upload.sh)) |

### DeJager FASTQs

| | |
|---|---|
| **Description** | ~50 multiplexed libraries, ~9.8 TB |
| **Primary location** | Synapse project [`syn21438684`](https://www.synapse.org/#!Synapse:syn21438684) |
| **Requirements** | Synapse account + project access + personal access token |
| **Canonical path (Openmind)** | `$DEJAGER_FASTQS_DIR` = `Data/Transcriptomics/DeJager/FASTQs/` |
| **NAS backup** | `LabShared/mabdel03/ROSMAP/Data/Transcriptomics/DeJager/FASTQs/` |

**How to get DeJager FASTQs:**

| Scenario | Script |
|----------|--------|
| Download from Synapse | [`Preprocessing/DeJager/01_FASTQ_Download/Download_FASTQs.sh`](../Preprocessing/DeJager/01_FASTQ_Download/Download_FASTQs.sh) |
| NAS -> local | [`Tsai_Server/DeJager/download_from_nas.sh`](Transcriptomics/Tsai_Server/DeJager/download_from_nas.sh) |
| Engaging <-> Openmind | [`Engaging-Openmind_Transfer/DeJager/send_globus.sh`](Transcriptomics/Engaging-Openmind_Transfer/DeJager/send_globus.sh) / [`receive_globus.sh`](Transcriptomics/Engaging-Openmind_Transfer/DeJager/receive_globus.sh) |

**Synapse setup:**
```bash
source config/paths.sh && init_conda && conda activate "${SYNAPSE_ENV}"
synapse login --auth-token <YOUR_TOKEN>
sbatch Preprocessing/DeJager/01_FASTQ_Download/Download_FASTQs.sh
```

---

## Transcriptomics -- Preprocessing Outputs

### CellRanger Output

Generated on Openmind from FASTQs using Cell Ranger 8.0.0 (`$CELLRANGER_PATH`).

| | Tsai | DeJager |
|---|---|---|
| **Path variable** | `$TSAI_CELLRANGER` | `$DEJAGER_CELLRANGER` |
| **Location** | `Data/Transcriptomics/Tsai/Cellranger_Output/` | `Data/Transcriptomics/DeJager/Cellranger_Output/` |
| **Samples** | ~480 | ~50 |
| **Key files per sample** | `raw_feature_bc_matrix.h5`, `filtered_feature_bc_matrix.h5`, BAM | Same |
| **NAS backup** | `LabShared/.../Tsai/Cellranger_Output/` | `LabShared/.../DeJager/Cellranger_Output/` |
| **Preprocessing scripts** | [`Preprocessing/Tsai/02_Cellranger_Counts/`](../Preprocessing/Tsai/02_Cellranger_Counts/) | [`Preprocessing/DeJager/02_Cellranger_Counts/`](../Preprocessing/DeJager/02_Cellranger_Counts/) |

Transfer scripts: same as FASTQs (NAS upload/download and Globus send/receive in the corresponding `Tsai/` or `DeJager/` subdirectories).

### CellBender Output

Ambient RNA removal applied to CellRanger output.

| | Tsai | DeJager |
|---|---|---|
| **Path variable** | `$TSAI_CELLBENDER` | `$DEJAGER_CELLBENDER` |
| **Location** | `Data/Transcriptomics/Tsai/Cellbender_Output/` | `Data/Transcriptomics/DeJager/Cellbender_Output/` |
| **Samples** | 478 | ~127 |
| **Size** | ~862 GB | ~525 GB |
| **Key file per sample** | `processed_feature_bc_matrix_filtered.h5` | Same |
| **Preprocessing scripts** | [`Preprocessing/Tsai/03_Cellbender/`](../Preprocessing/Tsai/03_Cellbender/) | [`Preprocessing/DeJager/03_Cellbender/`](../Preprocessing/DeJager/03_Cellbender/) |

**NAS backup notes:**
- Tsai CellBender is also at the legacy NAS path: `LabShared/mabdel03/ROSMAP/RNAseq/Tsai_Sequencing/Cellbender_Outputs` (478 folders, uploaded March 2026)
- DeJager CellBender primary copy is on permanent storage: `$DATA_ROOT/Data/DeJager/Preprocessed_Counts/`

---

## Transcriptomics -- Processing Outputs

The processing pipeline has 3 stages. Each stage's output is the input to the next.

Processing outputs are backed up on the NAS and can be downloaded instead of regenerated:

```bash
# Download all processing outputs (QC-filtered, doublet-removed, annotated)
cd Data_Access/Transcriptomics/Tsai_Server/Tsai
./download_processing_outputs_from_nas.sh

# Download only the annotated h5ad (Stage 3 output)
./download_processing_outputs_from_nas.sh 03_Integrated

# Upload processing outputs to NAS for backup
./upload_processing_outputs_to_nas.sh upload
```

To regenerate from scratch instead:
```bash
cd Processing/{Tsai,DeJager}/Pipeline/
./submit_pipeline.sh <stage>     # e.g., ./submit_pipeline.sh 1
./submit_pipeline.sh all         # all 3 stages with SLURM dependency chaining
```

### Stage 1: QC-Filtered

Per-sample quality control filtering (mitochondrial %, gene counts, total counts).

| | Tsai | DeJager |
|---|---|---|
| **Path variable** | `$TSAI_QC_FILTERED` | `$DEJAGER_QC_FILTERED` |
| **Output per sample** | `{sample}_qc.h5ad` | Same |
| **Summary** | `qc_summary.csv` | Same |
| **Script** | [`Processing/Tsai/Pipeline/01_qc_filter.py`](../Processing/Tsai/Pipeline/01_qc_filter.py) | [`Processing/DeJager/Pipeline/01_qc_filter.py`](../Processing/DeJager/Pipeline/01_qc_filter.py) |
| **SLURM** | Array job, 4 cores, 32 GB, ~1h/sample | Same |
| **Conda env** | `$QC_ENV` | Same |

### Stage 2: Doublet-Removed

Computational doublet detection and removal using scDblFinder.

| | Tsai | DeJager |
|---|---|---|
| **Path variable** | `$TSAI_DOUBLET_REMOVED` | `$DEJAGER_DOUBLET_REMOVED` |
| **Output per sample** | `{sample}_singlets.h5ad` | Same |
| **Summary** | `doublet_summary.csv` | Same |
| **Script** | [`Processing/Tsai/Pipeline/02_doublet_removal.Rscript`](../Processing/Tsai/Pipeline/02_doublet_removal.Rscript) | [`Processing/DeJager/Pipeline/02_doublet_removal.Rscript`](../Processing/DeJager/Pipeline/02_doublet_removal.Rscript) |
| **SLURM** | Array job, 4 cores, 32 GB, ~1h/sample | Same |
| **Conda env** | `$SINGLECELL_ENV` | Same |

### Stage 3: Integrated & Cell-Type-Annotated

Merges all samples, runs Harmony batch correction, Leiden clustering, UMAP, and cell type annotation via ORA with Mohammadi 2020 PFC markers.

| | Tsai | DeJager |
|---|---|---|
| **Path variable** | `$TSAI_INTEGRATED` | `$DEJAGER_INTEGRATED` |
| **Main output** | `tsai_annotated.h5ad` (~83 GB) | `dejager_annotated.h5ad` |
| **Script** | [`Processing/Tsai/Pipeline/03_integration_annotation.py`](../Processing/Tsai/Pipeline/03_integration_annotation.py) | [`Processing/DeJager/Pipeline/03_integration_annotation.py`](../Processing/DeJager/Pipeline/03_integration_annotation.py) |
| **SLURM** | Single job, 32 cores, 500 GB, 4-8h | Same |
| **Conda env** | `$BATCHCORR_ENV` | Same |

**Key fields in the annotated object (`obs`):**

| Field | Description |
|-------|-------------|
| `cell_type` | Top-1 cell type annotation per Leiden 0.5 cluster |
| `leiden_res0_2`, `leiden_res0_5`, `leiden_res1_0` | Cluster assignments at 3 resolutions |
| `patient_id` / `projid` | Patient identifier |
| `n_genes_by_counts`, `total_counts`, `pct_counts_mt` | QC metrics |

**Additional Stage 3 outputs:**
- `cluster_annotation_rankings.csv` -- ORA enrichment scores for all cell types per cluster
- `cluster_annotation_top3.csv` -- top 3 cell type annotations per cluster
- `figures/` -- UMAPs, elbow plots, batch correction comparisons, cell type proportions, marker dotplots

**Marker reference:** [`Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds`](../Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds)

---

## Phenotypes / Clinical Data

Small CSVs tracked in the repo at `Data/Phenotypes/` (`$PHENOTYPE_DIR`).

| File | Description | Rows |
|------|-------------|------|
| `ROSMAP_clinical.csv` | Core ROSMAP clinical variables (sex, education, APOE, cogdx, Braak, CERAD, PMI) | 3,584 |
| `dataset_652_basic_12-23-2021.csv` | Comprehensive phenotype extract (cognitive, neuropath, biomarkers, lifestyle) | 3,681 |
| `TSAI_DEJAGER_all_patients_wACEscores.csv` | Tsai + DeJager patients with ACE scores | 296 |
| `DeJager_ID_Map.csv` | Maps `projid` to `individualID` for DeJager samples | ~20 |

**NAS backup:** `LabShared/mabdel03/ROSMAP/Data/Phenotypes/`

**Transfer scripts:**
- Upload to NAS: [`Phenotypes/Tsai_Server/upload_to_nas.sh`](Phenotypes/Tsai_Server/upload_to_nas.sh)
- Download from NAS: [`Phenotypes/Tsai_Server/download_from_nas.sh`](Phenotypes/Tsai_Server/download_from_nas.sh)
- Globus send/receive: [`Phenotypes/Engaging-Openmind_Transfer/`](Phenotypes/Engaging-Openmind_Transfer/)

---

## Reference Files

| Resource | Path Variable | Location |
|----------|--------------|----------|
| Cell Ranger 8.0.0 | `$CELLRANGER_PATH` | `/om2/user/mabdel03/apps/yard/cellranger-8.0.0` |
| GRCh38 reference | `$CELLRANGER_REF` | `/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A` |
| PFC marker genes | `$TSAI_MARKERS_RDS` | `Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds` |
| Tsai patient metadata | `$TSAI_METADATA_CSV` | `Preprocessing/Tsai/02_Cellranger_Counts/Tracking/patient_metadata.csv` |
| DeJager barcode-to-patient map | `$DEJAGER_PATIENT_MAP` | `/om/scratch/Mon/shared_folder/WGS/cell_to_patient_assignmentsFinal0.csv` |
| DeJager patient ID overrides | `$DEJAGER_PATIENT_ID_OVERRIDES` | `Processing/DeJager/Pipeline/Resources/patient_id_overrides.json` |
| All path variables | -- | [`config/paths.sh`](../config/paths.sh) |

---

## Transfer Methods Reference

### Tsai Server (TsaiLabNAS) -- SFTP

- **Protocol**: SFTP with `sshpass` for non-interactive authentication
- **Host**: `tsailabnas.mit.edu` (Synology NAS, DSM 7.3.2)
- **Authentication**: Password file at `~/.smb_tsailabnas` (SMB format: `password = <pass>`)
- **NAS base path**: `/LabShared/mabdel03/ROSMAP/`
- **Parallelism**: SLURM array jobs with SFTP sharding for large transfers
- **Constraints**: No SSH shell access on NAS; rsync daemon auth disabled; smbclient not available on compute nodes. SFTP is the only working protocol.

### Engaging-Openmind Transfer -- Globus

- **Protocol**: Globus (recommended for large transfers)
- **GUI**: https://app.globus.org/
- **CLI**: `globus transfer` command (installed in conda environments)
- **Endpoint IDs**:
  - Openmind: `cbc6f8da-d37e-11eb-bde9-5111456017d9`
  - Engaging: `c52fcff2-761c-11eb-8cfc-cd623f92e1c0`
- **Conda environments**:
  - Openmind: `/om2/user/mabdel03/conda_envs/globus_env`
  - Engaging: `/home/mabdel03/conda_envs/globus_env`

### Prerequisites

#### For Tsai Server transfers:
1. Password file at `~/.smb_tsailabnas` with format:
   ```
   username = mabdel03
   password = YOUR_PASSWORD
   ```
2. `chmod 600 ~/.smb_tsailabnas`
3. `sshpass` available (installed system-wide on OpenMind)

#### For Globus transfers:
1. Activate the Globus conda environment
2. Run `globus login` (one-time setup, opens browser)
3. Verify with `globus whoami`

### Data Layout on NAS

NAS paths mirror the local `Data/` structure:
```
LabShared/mabdel03/ROSMAP/Data/
├── Transcriptomics/
│   ├── Tsai/
│   │   ├── FASTQs/
│   │   ├── Cellranger_Output/
│   │   ├── Cellbender_Output/
│   │   └── Processing_Outputs/
│   │       ├── 01_QC_Filtered/
│   │       ├── 02_Doublet_Removed/
│   │       └── 03_Integrated/
│   └── DeJager/
│       ├── FASTQs/
│       ├── Cellranger_Output/
│       ├── Cellbender_Output/
│       └── Processing_Outputs/
│           ├── 01_QC_Filtered/
│           ├── 02_Doublet_Removed/
│           └── 03_Integrated/
└── Phenotypes/
```

Note: Tsai CellBender data was previously uploaded to
`LabShared/mabdel03/ROSMAP/RNAseq/Tsai_Sequencing/Cellbender_Outputs` (478 folders).

### Transfer Scripts Directory Structure

```
Data_Access/
├── README.md                                    (this file)
├── Transcriptomics/
│   ├── Tsai_Server/                             SFTP transfers to/from TsaiLabNAS
│   │   ├── README.md
│   │   ├── Tsai/
│   │   │   ├── upload_to_nas.sh
│   │   │   ├── download_from_nas.sh
│   │   │   ├── upload_to_nas.env
│   │   │   ├── slurm_upload.sh
│   │   │   ├── upload_processing_outputs_to_nas.sh
│   │   │   ├── download_processing_outputs_from_nas.sh
│   │   │   └── upload_processing_outputs_to_nas.env
│   │   └── DeJager/
│   │       ├── upload_to_nas.sh
│   │       ├── download_from_nas.sh
│   │       ├── upload_to_nas.env
│   │       ├── slurm_upload.sh
│   │       ├── upload_processing_outputs_to_nas.sh
│   │       ├── download_processing_outputs_from_nas.sh
│   │       └── upload_processing_outputs_to_nas.env
│   └── Engaging-Openmind_Transfer/              Globus transfers between HPC clusters
│       ├── README.md
│       ├── Tsai/
│       │   ├── generate_batch.sh
│       │   ├── send_globus.sh
│       │   └── receive_globus.sh
│       └── DeJager/
│           ├── generate_batch.sh
│           ├── send_globus.sh
│           └── receive_globus.sh
└── Phenotypes/
    ├── Tsai_Server/
    │   ├── upload_to_nas.sh
    │   └── download_from_nas.sh
    └── Engaging-Openmind_Transfer/
        ├── send_globus.sh
        └── receive_globus.sh
```
