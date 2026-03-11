# Scientific Background

This document explains the scientific context behind the ROSMAP snRNA-seq pipeline for readers who may be unfamiliar with the dataset, the technology, or the rationale behind each processing step.

## What is ROSMAP?

ROSMAP refers to two complementary longitudinal cohort studies of aging and Alzheimer's disease (AD):

- **Religious Orders Study (ROS)** — started in 1994, enrolling older Catholic nuns, priests, and brothers
- **Memory and Aging Project (MAP)** — started in 1997, enrolling older adults from the greater Chicago area

Both studies are conducted by Rush University Medical Center. Participants undergo annual cognitive and clinical evaluations, and all have agreed to organ donation at the time of death. This provides matched antemortem clinical data (cognitive trajectories, clinical diagnoses, lifestyle factors) with postmortem brain tissue, making ROSMAP one of the most deeply phenotyped brain tissue repositories available for AD research.

## Why the Dorsolateral Prefrontal Cortex (DLPFC)?

This pipeline processes tissue from the **dorsolateral prefrontal cortex (DLPFC)**, corresponding roughly to Brodmann areas 9 and 46. The DLPFC is one of the most commonly profiled brain regions in AD research because:

- It is involved in higher-order cognitive functions (working memory, executive function) that decline in AD
- It is affected across a broad range of disease severity, from mild cognitive impairment to late-stage dementia
- Unlike the entorhinal cortex or hippocampus (which show pathology earliest), the DLPFC captures a wider spectrum of disease states within a single cohort

## Why Single-Nucleus RNA Sequencing (snRNA-seq)?

Traditional single-cell RNA sequencing (scRNA-seq) requires fresh, dissociated tissue, which is impractical for postmortem brain samples. **Single-nucleus RNA sequencing** isolates individual nuclei from frozen tissue, enabling transcriptomic profiling at single-cell resolution from archived biobank samples. This is critical for ROSMAP, where tissue is collected at autopsy and stored frozen.

snRNA-seq captures the transcriptional state of each cell type in the brain (neurons, astrocytes, microglia, oligodendrocytes, etc.), enabling:

- Cell type-specific differential expression analysis
- Identification of disease-associated cell states
- Characterization of cell type proportions across disease conditions

## Why `--include-introns true`?

Nuclear RNA has a different composition than cytoplasmic RNA. Because snRNA-seq captures RNA from the nucleus (before most splicing is complete), a substantial fraction of reads map to **intronic regions** rather than exons. Including intronic reads with `--include-introns true` in Cell Ranger significantly increases the number of detected genes and UMIs per nucleus, improving statistical power for downstream analysis. Without this flag, a large portion of genuinely informative reads would be discarded.

## Why CellBender (Ambient RNA Removal)?

During nucleus isolation, some RNA leaks out of damaged or lysed cells into the surrounding solution. This **ambient RNA** is then captured in droplets alongside intact nuclei, contaminating the expression profile of every cell with a background signal that reflects the overall tissue composition rather than individual cell identity.

CellBender uses a deep generative model to distinguish true cell-associated RNA from this ambient background. Without this correction:

- Cell type markers can "bleed" across cell types (e.g., neuronal markers appearing in non-neuronal cells)
- Differential expression results can be confounded by ambient contamination levels rather than true biological differences
- Clustering and cell type annotation become less accurate

The pipeline uses `--fpr 0` (or `0.01`) for a stringent false positive rate, aggressively removing ambient signal.

## Why MAD-Based QC Filtering?

After ambient RNA removal, individual cells are filtered to remove low-quality observations. The pipeline uses **Median Absolute Deviation (MAD)** as a robust measure of variability:

- **Why MAD instead of standard deviation?** MAD is less sensitive to extreme outliers, making it more appropriate for distributions with long tails (common in single-cell data).
- **What is filtered?** Cells with metrics that fall beyond N MADs from the median are flagged as outliers:
  - `log1p_total_counts` (4 MADs) — removes cells with abnormally high or low total RNA counts
  - `log1p_n_genes_by_counts` (4 MADs) — removes cells expressing too many or too few genes
  - `pct_counts_in_top_20_genes` (4 MADs) — removes cells dominated by a handful of genes
  - `pct_counts_mt` (3 MADs or >7.5%) — removes cells with high mitochondrial content, a marker of cell damage

These thresholds balance removing damaged/dying cells against retaining genuine biological variation. The 4-MAD threshold for counts and genes is relatively permissive; the stricter 3-MAD + 7.5% cap for mitochondrial percentage reflects the stronger association between mitochondrial content and cell death.

## Why Harmony for Batch Correction?

The Tsai dataset spans 16 processing batches (30 patients each). Even with identical protocols, technical differences between batches (reagent lots, sequencing runs, handling variation) introduce systematic variation that can confound biological signal.

**Harmony** corrects for batch effects by iteratively adjusting the PCA embedding so that cells from different batches overlap in reduced-dimensional space. Critically:

- It operates on the **PCA embedding** only, leaving the raw count matrix untouched
- This preserves the original expression values for downstream differential expression analysis
- The correction uses the patient ID as the covariate (`projid` for Tsai, `patient_id` for DeJager), treating each patient as a "batch" to remove per-patient technical variation while preserving condition-level biological signal

## The Two Datasets

### DeJager Dataset

- **Source**: Downloaded from Synapse (project syn21438684)
- **Key characteristic**: Libraries are **multiplexed** — each sequencing library contains cells from multiple patients pooled together
- **Extra step required**: Genotype-based demultiplexing (Demuxlet/Freemuxlet) using whole-genome sequencing (WGS) data to assign each cell barcode to its patient of origin
- **Pipeline**: Download FASTQs from Synapse, Cell Ranger alignment, CellBender, then Demuxlet

### Tsai Dataset

- **Source**: Located on the MIT Engaging cluster
- **Key characteristic**: Patient assignments are **known** from sequencing metadata (one patient per library)
- **No demultiplexing needed**: Cell-to-patient assignment comes directly from the library metadata
- **Scale**: 480 patients, 5,197 FASTQ files (~9 TB), processed in 16 batches of 30 patients
- **Cohorts**: Patients are drawn from three study arms:
  - **ACE** — Adverse Childhood Experiences
  - **Resilient** — Cognitive Resilience (maintained cognitive function despite AD pathology)
  - **SocIsl** — Social Isolation

Both datasets are processed through the same core pipeline (Cell Ranger, CellBender, QC, doublet removal, integration, annotation) to ensure comparability.

## What the Pipeline Produces

The final output is an **annotated AnnData object** (`tsai_annotated.h5ad`) containing:

- Single-cell gene expression data for all patients, filtered and quality-controlled
- Cell type annotations (assigned via over-representation analysis using Mohammadi 2020 brain cell type markers)
- Batch-corrected embeddings (Harmony-adjusted PCA and UMAP)
- Patient and batch metadata for each cell

This object is the starting point for **downstream analysis**:

- **Differential expression (DEG)**: Compare gene expression between conditions (AD vs. control, resilient vs. non-resilient) within specific cell types using pseudobulk methods
- **SCENIC**: Infer gene regulatory networks and identify active transcription factor regulons per cell type
- **Transcription factor analysis**: Characterize TF activity differences across disease states
