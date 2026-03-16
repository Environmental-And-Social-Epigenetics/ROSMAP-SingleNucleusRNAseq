# SocIsl — Social Isolation

Analysis of the impact of social isolation on gene expression in Alzheimer's disease brain tissue.

## Phenotype Definition

Social isolation is measured using ROSMAP variables:
- **Social network size** (`soc_net_bl`, `soc_net_lv`)
- **Social isolation score** (`social_isolation_avg`, `social_isolation_lv`)

## Patient Selection

Phenotype data: `Data/Phenotypes/dataset_652_basic_12-23-2021.csv`

| Group | Criteria | Description |
|-------|----------|-------------|
| Socially isolated | High `social_isolation_avg` | Patients with elevated social isolation scores |
| Non-isolated | Low `social_isolation_avg` | Patients with typical social engagement |

## Key Comparisons

- Isolated vs non-isolated, stratified by cell type and sex
- Continuous association of social isolation score with gene expression
- Interaction with AD pathology

## Status

This analysis was originally conducted on Openmind (2024-2025) using both Tsai and DeJager datasets. The scripts and small result files have been migrated into this repository structure. Large data outputs (.h5ad, .rds pseudobulk matrices, .tsv expression matrices, feather ranking files) are backed up on the Tsai Lab NAS.

## Directory Structure

```
SocIsl/
├── DEG/                    Differential expression (limma/DESeq2)
│   ├── DeJager/            Scripts + limma results (12 cell type x sex CSVs)
│   └── Tsai/               Scripts for Tsai DEG analysis
├── SCENIC/                 Regulatory network inference (pySCENIC)
│   ├── DeJager/            SCENIC scripts for DeJager
│   └── Tsai/               SCENIC scripts + regulon/AUCell result CSVs
├── TF/                     Transcription factor / COMPASS metabolic analysis
│   ├── DeJager/            COMPASS run scripts (per cell type) + cell ranger scripts
│   └── Tsai/               COMPASS run scripts (per cell type x sex) + TF activity CSVs
├── GSEA/                   Gene set enrichment / pathway analysis (WebGestaltR)
│   ├── DeJager/            GSEA scripts + pathway result CSVs
│   └── Tsai/               Tsai-specific GSEA scripts
├── _data_prep/             Legacy data preparation and processing scripts
│   ├── DeJager/            DeJager preprocessing, library matching, data download
│   └── Tsai/               Tsai preprocessing, matrix extraction, variable setup
└── figures/                Output plots (volcano, heatmap, RRHO, UMAP, COMPASS)
```

## Analysis Types

### DEG (Differential Expression)
- **Method**: limma (primary), DESeq2 (secondary), Nebula (mixed model)
- **Design**: Social isolation score as continuous covariate, stratified by sex and cell type
- **Scripts**: `socIslDeg.Rscript` (main DEG), `dejagerNew.Rscript` (DeJager replication)
- **Results**: `limma_SIA_results{Sex}{CellType}.csv` (12 files per dataset)

### SCENIC (Single-Cell Regulatory Networks)
- **Method**: pySCENIC with hg38 motif rankings
- **Reference files needed** (not in repo, ~3.5GB total):
  - `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
- **Results**: `{sex}_regulonsFULL_Tsai{CellType}.csv`, `{sex}_auc_mtxFULL_Tsai{CellType}.csv`

### TF / COMPASS (Metabolic Flux Analysis)
- **Method**: COMPASS (metabolic modeling) with CPLEX solver
- **Note**: Requires IBM CPLEX license (academic license available)
- **Design**: Per-cell-type, per-sex COMPASS runs on both datasets
- **Results**: p-values (`{sex}_pVals{CellType}.csv`), metadata (`{sex}_metaDF{CellType}.csv`)
- **Engaging location**: DeJager COMPASS outputs at `/home/nkhera/orcd/pool/Subfolder/` (507GB)

### GSEA (Gene Set Enrichment)
- **Method**: WebGestaltR (ORA)
- **Pathways**: KEGG, Reactome, Panther, Wikipathway, GO (BP, CC, MF)
- **Results**: `plotGsea{Sex}{Pathway}.csv`

## Reproducing This Analysis on Engaging

### Prerequisites

1. Complete the Processing pipeline first (`Processing/Tsai/Pipeline/submit_pipeline.sh all`)
   to generate `tsai_annotated.h5ad`
2. Install Analysis conda environments: `bash setup/install_envs.sh --analysis`
3. Ensure phenotype data is in `Data/Phenotypes/dataset_652_basic_12-23-2021.csv`

### Step-by-Step

```bash
source config/paths.sh

# --- DEG Analysis (simplest starting point) ---
init_conda && conda activate "${CONDA_ENV_BASE}/deg_analysis"
cd Analysis/SocIsl/DEG/Tsai
# Edit socIslDegT.sh: update SLURM headers for Engaging partitions
sbatch socIslDegT.sh
# Output: limma_SIA_results{Sex}{CellType}.csv

# --- GSEA (requires DEG results) ---
conda activate "${CONDA_ENV_BASE}/gsea_analysis"
cd Analysis/SocIsl/GSEA/DeJager
# Edit gsea.sh: update SLURM headers
sbatch gsea.sh
# Output: plotGsea{Sex}{Pathway}.csv

# --- SCENIC (requires annotated h5ad, long-running) ---
conda activate "${CONDA_ENV_BASE}/scenic_analysis"
# Download reference files first (~3.5 GB):
#   wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#   wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#   wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
cd Analysis/SocIsl/SCENIC/Tsai
sbatch tsaiAdataScenic.sh

# --- COMPASS / TF (requires CPLEX license) ---
conda activate "${CONDA_ENV_BASE}/compass_analysis"
# Install CPLEX: https://www.ibm.com/academic/ (free academic license)
cd Analysis/SocIsl/TF/Tsai
sbatch compassRunAstTsai.sh   # one cell type at a time
```

### About `_data_prep/`

The `_data_prep/` scripts are **NOT needed** if you have the annotated h5ad from
the Processing pipeline. They were the original ad-hoc preprocessing workflow
used before the standardized Processing pipeline existed. They are preserved
for historical reference only.

## Updating Hardcoded Paths

Legacy scripts contain hardcoded Openmind paths. Before running any script,
update these path patterns:

| Old Path (Openmind) | Replacement (Engaging) |
|---------------------|------------------------|
| `/om2/user/mabdel03/anaconda/etc/profile.d/conda.sh` | `source config/paths.sh && init_conda` (or your `$CONDA_INIT_SCRIPT`) |
| `/om2/user/mabdel03/conda_envs/<env>` | `${CONDA_ENV_BASE}/<env>` |
| `/om/scratch/Mon/mabdel03/SocialIsolation/` | Your working directory on Engaging |
| `/om/scratch/Sun/mabdel03/SocialIsolation/` | Your working directory on Engaging |
| `/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/` | `${DATA_ROOT}/` |
| `/om/scratch/Mon/mabdel03/SocialIsolation/opt/ibm/ILOG/CPLEX_Studio2211` | Your CPLEX install path on Engaging |

To do a bulk find-and-replace preview across all scripts:

```bash
grep -rn "/om2/\|/om/scratch\|/net/vast" Analysis/SocIsl/ --include="*.sh" --include="*.py" --include="*.Rscript"
```

## Important Notes

- Large data outputs (pseudobulk .tsv matrices, .rds objects, .h5ad files, COMPASS output directories ~507GB) are **not included** in this repo — they are backed up on the Tsai Lab NAS.
- DeJager COMPASS outputs are on Engaging at `/home/nkhera/orcd/pool/Subfolder/` (507GB).
