# Resources

Reference data files used by the processing pipeline.

## Files

### `Brain_Human_PFC_Markers_Mohammadi2020.rds`

Cell type marker gene set from Mohammadi et al. (2020), containing marker genes for human prefrontal cortex (PFC) cell types. This RDS file is loaded by Stage 3 of the pipeline (`03_integration_annotation.py`) for over-representation analysis (ORA)-based cell type annotation.

The marker set includes genes for major brain cell types (excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, microglia, OPCs, etc.) and is used to assign cluster identities after Harmony integration.

## Reference

Mohammadi, S., Davila-Velderrain, J., & Bhatt, P. (2020). A multiresolution framework to characterize single-cell state landscapes. *Nature Communications*, 11, 5399.
