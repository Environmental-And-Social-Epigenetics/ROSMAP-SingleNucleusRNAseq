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

## Important Notes

- These scripts contain **hardcoded paths from the Openmind cluster** and are preserved for reference. They will need path updates to run on Engaging.
- Large data outputs (pseudobulk .tsv matrices, .rds objects, .h5ad files, COMPASS output directories) are **not included** in this repo — they are backed up on the Tsai Lab NAS.
- The `_data_prep/` scripts show the original data preparation workflow but are superseded by the current Processing pipeline.
