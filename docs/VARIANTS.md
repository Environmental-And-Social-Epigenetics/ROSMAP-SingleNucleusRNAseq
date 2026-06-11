# Processing Variants

Variant IDs are defined in `config/variants.yaml` and are the only supported names for Stage 3 integration variants.

## Tsai

| Variant | Role | Harmony key | Notes |
|---------|------|-------------|-------|
| `derived_batch` | canonical | `derived_batch` | Flowcell-derived batch correction. |
| `projid` | sensitivity | `projid` | Earlier patient-level correction. |
| `no_harmony` | sensitivity | none | PCA baseline with Harmony skipped. |

## DeJager

| Variant | Role | Harmony key | Notes |
|---------|------|-------------|-------|
| `library_id` | canonical | `library_id` | Per-library correction. |
| `patient_id` | sensitivity | `patient_id` | Patient-level correction. |
| `pool_batch` | sensitivity | `pool_batch` | B-number pool correction. |
| `derived_batch` | sensitivity | `derived_batch` | Flowcell-derived correction. |
| `sequencing_date` | sensitivity | `sequencing_date` | YYMMDD library-prefix correction. |
| `no_harmony` | sensitivity | none | PCA baseline with Harmony skipped. |

Outputs are written to `Processing_Outputs/03_Integrated/{variant_id}/`. Compatibility variables in `config/paths.sh` point to these canonical variant paths.

