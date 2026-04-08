# ACE TF/COMPASS — Metabolic Flux Analysis

COMPASS metabolic flux analysis for the ACE phenotype, adapted from the SocIsl
implementation in `Analysis/SocIsl/TF/`.

## Prerequisites

**IBM CPLEX 2211** (academic license required):
1. Register at https://www.ibm.com/academic/
2. Download CPLEX Studio 2211
3. Set in `config/paths.local.sh`:
   ```bash
   export CPLEX_DIR="/path/to/ibm/ILOG/CPLEX_Studio2211"
   ```
4. Install Python bindings: `pip install cplex`

## Method

1. Generate pseudobulk TSV expression matrices per cell type (from DEG splits)
2. Run COMPASS metabolic flux prediction using CPLEX solver
3. Convert reaction penalties to activity scores
4. Mann-Whitney U test between ACE phenotype groups
5. FDR correction and effect size calculation

## Running

```bash
source config/paths.sh
cd Analysis/ACE/TF/Tsai
sbatch compassRun_{celltype}.sh
```

After COMPASS runs complete, run post-processing:
```bash
python compass_analysis.py
```

## Outputs

Under `${ACE_OUTPUT_ROOT}/TF/Tsai/`:
- `{celltype}Compass/` — raw COMPASS outputs per cell type
- Reaction activity matrices, meta-reaction assignments
- Statistical results with p-values, effect sizes, adjusted p-values

## Subdirectories

- `DeJager/` — TF analysis on DeJager dataset
- `Tsai/` — TF analysis on Tsai dataset
