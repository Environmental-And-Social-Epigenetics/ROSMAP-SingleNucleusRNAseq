# Phenotypes

Clinical phenotype data from the ROSMAP dataset. These files map patient IDs (`projid`) to clinical, demographic, neuropathological, and behavioral phenotypes.

## Files

### ROSMAP_clinical.csv

| Property | Value |
|----------|-------|
| Rows | 3,584 |
| Size | 257 KB |
| Source | ROSMAP clinical core dataset |

Core clinical data for all ROSMAP participants.

**Key columns**: `projid`, `Study`, `msex`, `educ`, `race`, `apoe_genotype`, `age_death`, `cogdx` (cognitive diagnosis), `braaksc` (Braak stage), `ceradsc` (CERAD score), `pmi`, `individualID`

**Used by**: All phenotype analyses (baseline demographics and AD diagnosis)

---

### dataset_652_basic_12-23-2021.csv

| Property | Value |
|----------|-------|
| Rows | 3,681 |
| Size | 1.8 MB |
| Source | ROSMAP dataset 652 (comprehensive phenotype extract) |

Comprehensive phenotype file with ACE scores, cognitive measures, neuropathology, biomarkers, personality traits, and lifestyle variables.

**Key columns**:
- **ACE**: `emotional_neglect`, `family_pro_sep`, `financial_need`, `parental_intimidation`, `parental_violence`, `tot_adverse_exp`
- **Cognitive**: `cogdx`, `cogn_global_lv`, `cogn_ep_lv`, `cogn_po_lv`, `cogn_ps_lv`, `cogn_se_lv`, `cogn_wo_lv`
- **Neuropathology**: `braaksc`, `ceradsc`, `amyloid`, `plaq_d`, `plaq_n`, `nft`, `tangles`
- **Social**: `soc_net_bl`, `social_isolation_avg`, `social_isolation_lv`
- **Personality**: `agreeableness`, `conscientiousness`, `extraversion_6`, `neuroticism_12`, `openness`
- **Biomarkers**: inflammatory (IL-6, TNFa, VCAM), metabolic (glucose, cholesterol, HbA1c)

**Used by**: ACE analysis, SocIsl analysis, Resilient analysis

---

### TSAI_DEJAGER_all_patients_wACEscores.csv

| Property | Value |
|----------|-------|
| Rows | 296 |
| Size | 656 KB |
| Source | Filtered to Tsai + DeJager patients with ACE data |

Subset of patients present in both the snRNA-seq datasets (Tsai and DeJager) who have ACE scores. Contains all columns from the comprehensive dataset plus additional biomarker timepoints (baseline, last visit, average).

**Key columns**: Same as `dataset_652_basic_12-23-2021.csv` plus longitudinal biomarker measures (`_bl`, `_lv`, `_avg` suffixes)

**Note**: `projid` is the **last** column in this file. Row index is the first column.

**Used by**: ACE analysis (primary phenotype file for ACE comparisons)

---

### DeJager_ID_Map.csv

| Property | Value |
|----------|-------|
| Rows | ~20 |
| Size | 1.5 KB |
| Source | Manual mapping |

Maps `projid` to `individualID` for the DeJager dataset. Required because DeJager uses `individualID` in some metadata files while the pipeline uses `projid`.

**Columns**: row index, `projid`, `individualID`

## Config Variables

These files are referenced in `config/paths.sh`:

```bash
PHENOTYPE_DIR="${REPO_ROOT}/Data/Phenotypes"
ROSMAP_CLINICAL_CSV="${PHENOTYPE_DIR}/ROSMAP_clinical.csv"
ACE_SCORES_CSV="${PHENOTYPE_DIR}/TSAI_DEJAGER_all_patients_wACEscores.csv"
DEJAGER_ID_MAP="${PHENOTYPE_DIR}/DeJager_ID_Map.csv"
```

## Provenance

All files originate from the ROSMAP data portal (Synapse) and were extracted/filtered for this project. The comprehensive dataset (`dataset_652_basic_12-23-2021.csv`) was downloaded on 2021-12-23.
