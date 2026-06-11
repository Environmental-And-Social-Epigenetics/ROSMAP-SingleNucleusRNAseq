# ACE Co-Expression Network Analysis (hdWGCNA)

Identifies co-expressed gene modules at single-cell resolution using hdWGCNA,
then tests module-trait associations with ACE phenotypes.

## Method

hdWGCNA extends WGCNA for single-cell data via metacell construction. For each
cell type, cells are aggregated into metacells (k=25), WGCNA identifies
co-expressed gene modules, and module eigengenes are correlated with ACE
phenotypes at the patient level.

## Cohort Directories

| Cohort | Entry Point |
|--------|-------------|
| Tsai | `Tsai/aceWgcnaT.sh` |
| DeJager | `DeJager/aceWgcnaDJ.sh` |

## What This Analysis Adds

With 2,200 DEGs, the narrative is currently gene-level. hdWGCNA groups these
into functional co-expression modules (e.g., "synaptic signaling module,"
"chromatin remodeling module"), enabling statements like "a complement/phagocytosis
module is suppressed in microglia" rather than "327 genes are downregulated."

Modules also enable cross-cell-type comparison: shared modules between
microglia and PV+ basket cells would indicate coordinated disruption.
