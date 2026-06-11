# ACE

ACE analysis studies adverse childhood experience phenotypes in the ROSMAP
single-nucleus datasets.

## Official ACE Phenotype Models

These are the formalized primary models for both cohorts:

| Model | Variable | Definition |
|-------|----------|------------|
| Total adversity | `tot_adverse_exp` | sum of the five ACE component scores |
| Early household SES | `early_hh_ses` | tracked continuous SES measure from the phenotype table |
| Aggregate score | `ace_aggregate` | `zscore(tot_adverse_exp) - zscore(early_hh_ses)` (deprivation-aligned) |

ACE component columns are also available for exploratory composition models:

- `emotional_neglect`
- `family_pro_sep`
- `financial_need`
- `parental_intimidation`
- `parental_violence`

## Shared Input Contract

### Processed AnnData

- `cell_type`
- Tsai: `sample_id` or `projid`
- DeJager: `patient_id`

### Phenotype CSV

Default: `${ACE_SCORES_CSV}`

Required columns:

- `projid`
- `tot_adverse_exp`
- `early_hh_ses`
- `msex`
- `age_death`
- `pmi`
- `niareagansc`
- the five ACE component columns above

## Supported ACE Workflows

| Workflow | Tsai | DeJager | Status |
|----------|------|---------|--------|
| DEG | `DEG/Tsai/aceDegT.sh` | `DEG/DeJager/aceDegDJ.sh` | Complete |
| Cell-type proportion | `CellTypeProportion/Tsai/acePropT.sh` | `CellTypeProportion/DeJager/acePropDJ.sh` | Complete |
| GSEA | `GSEA/Tsai/aceGseaT.sh` | `GSEA/DeJager/aceGseaDJ.sh` | Ready |
| TF Activity | `TFActivity/Tsai/aceTfActT.sh` | `TFActivity/DeJager/aceTfActDJ.sh` | Ready |
| SCENIC | `SCENIC/Tsai/aceScenicT.sh` | `SCENIC/DeJager/aceScenicDJ.sh` | Ready |
| CellChat | `CellChat/Tsai/aceCellChatT.sh` | `CellChat/DeJager/aceCellChatDJ.sh` | **Blocked** (see below) |
| hdWGCNA | `hdWGCNA/Tsai/aceWgcnaT.sh` | `hdWGCNA/DeJager/aceWgcnaDJ.sh` | **Blocked** (see below) |
| Mic State | `MicState/Tsai/aceMicStateT.sh` | `MicState/DeJager/aceMicStateDJ.sh` | Ready |
| Epigenomic Integration | `EpigenomicIntegration/Tsai/aceEpiIntT.sh` | — | Ready (no DeJager scaffold) |

### Known issues

**CellChat (Tsai)** — `zellkonverter::readH5AD` on the 83 GB `tsai_annotated.h5ad`
times out. Two paths forward: (a) chunk the input into per-cell-type-group
h5ads upstream, (b) convert to RDS via a Python preprocessor and read RDS
directly. Conda env `cellchat_analysis` is built and CellChat 2.2 imports
cleanly; the blocker is purely the input format.

**hdWGCNA (Tsai)** — `zellkonverter::readH5AD` on the per-cell-type splits
fails with `Dimnames[[2]] (161176) is not equal to Dim[2] (36601)` — a
sparse-matrix orientation/encoding mismatch between anndata 0.12 and
zellkonverter 1.16. Same workaround options as CellChat. Conda env
`wgcna_analysis` is built and hdWGCNA 0.4.11 imports cleanly.

**MicState DEG-prep gotcha (resolved)** — phenotype loading now goes through
[_shared/load_ace_phenotype.py](_shared/load_ace_phenotype.py) so
`ace_aggregate` is computed on demand.

### Workflow Descriptions

| Workflow | What it does | Key output |
|----------|-------------|------------|
| **DEG** | Pseudobulk DESeq2 differential expression per cell type × sex | DEG lists, volcano plots |
| **Cell-type proportion** | sccomp compositional analysis | Proportion effect sizes |
| **GSEA** | Ranked gene set enrichment (8 pathway databases) | Enriched pathways per cell type |
| **TF Activity** | decoupler TF/pathway activity inference (DoRothEA, CollecTRI, PROGENy) | Differentially active TFs |
| **SCENIC** | pySCENIC gene regulatory networks (GRNBoost2 + motif validation) | Regulons, TF-target networks |
| **CellChat** | Cell-cell communication differential analysis | Altered L-R interactions |
| **hdWGCNA** | Co-expression network modules | Gene modules, trait correlations |
| **Mic State** | Microglial activation state scoring (homeostatic/DAM/inflammatory) | State scores, trajectory |
| **Epigenomic Integration** | Multi-modal integration with ROSMAP H3K9ac/methylation | Epigenetic mark changes at DEG loci |

Canonical integration choices:

- **Tsai**: `derived_batch`
- **DeJager**: `library_id`

Sensitivity integrations can still be passed explicitly to the launchers.

## Output Layout

All ACE outputs are written beneath:

```text
${ANALYSIS_OUTPUT_ROOT}/ACE/
  DEG/
    Tsai/
    DeJager/
  CellTypeProportion/
    Tsai/
    DeJager/
  GSEA/
    Tsai/
    DeJager/
  TFActivity/
    Tsai/
    DeJager/
  SCENIC/
    Tsai/
    DeJager/
  CellChat/
    Tsai/
    DeJager/
  hdWGCNA/
    Tsai/
    DeJager/
  MicState/
    Tsai/
    DeJager/
  EpigenomicIntegration/
    Tsai/
```

The repo directories store code only.

## Analysis Narrative

The analyses are organized in layers that build on each other:

```
Layer 1: WHAT changes
  ├── DEG: 2,200 genes, male-biased, downregulation-dominant
  ├── GSEA: Which pathways are disrupted in each cell type
  └── PROGENy: Which signaling pathways change activity

Layer 2: WHO controls it
  ├── TF Activity (decoupler): Which TFs are responsible
  ├── SCENIC: Full regulatory networks with motif validation
  └── Cross-cell-type TF convergence

Layer 3: HOW cells interact
  └── CellChat: Altered ligand-receptor signaling

Layer 4: ARCHITECTURE of disruption
  ├── hdWGCNA: Co-expression programs and modules
  ├── Mic State: DAM vs homeostatic shift
  └── Cell proportions: No cell loss (pure transcriptomic)

Layer 5: MECHANISM
  └── Epigenomic integration: H3K9ac/methylation evidence
```

## Validation Layers

- Smoke: fixture-based checks that exercise the ACE prep and analysis entrypoints
  without requiring the full processed cohort objects
- Full: production-readiness checks against the canonical annotated cohort inputs

Preflight entrypoints:

- `bash config/preflight.sh ace-tsai-smoke`
- `bash config/preflight.sh ace-dejager-smoke`
- `bash config/preflight.sh ace-tsai-full`
- `bash config/preflight.sh ace-dejager-full`
