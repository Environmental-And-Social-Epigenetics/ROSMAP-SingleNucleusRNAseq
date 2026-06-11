#!/bin/bash
# Smoke test for ACE GSEA pipeline (Tsai cohort)
#
# Creates a minimal mock DESeq2 .rda file and runs gsea_analysis.R in --smoke
# mode to verify the pipeline is functional end-to-end.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

SMOKE_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/Smoke/Tsai/GSEA"
FIXTURE_DIR="${SMOKE_ROOT}/fixtures/tot_adverse_exp"
RESULTS_DIR="${SMOKE_ROOT}/results"
mkdir -p "${FIXTURE_DIR}" "${RESULTS_DIR}"

echo "=== ACE GSEA Smoke Test ==="
echo "Fixture dir: ${FIXTURE_DIR}"
echo "Results dir: ${RESULTS_DIR}"
echo ""

# ── Activate environment ─────────────────────────────────────────────────────

set +u
activate_env "${GSEA_ANALYSIS_ENV}"
set -u

# ── Create mock .rda fixture ────────────────────────────────────────────────

FIXTURE_RDA="${FIXTURE_DIR}/deseqAnalysisACE_tot_adverse_exp_Exc_Female.rda"

cat > "${SMOKE_ROOT}/create_fixture.R" << 'REOF'
# Create a minimal DESeqResults-like object for smoke testing
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
output_path <- args[1]

set.seed(42)
n_genes <- 500

# Generate mock gene expression data (2 conditions, 3 replicates each)
count_matrix <- matrix(
  as.integer(rpois(n_genes * 6, lambda = 100)),
  nrow = n_genes, ncol = 6
)
rownames(count_matrix) <- paste0("GENE", seq_len(n_genes))
colnames(count_matrix) <- paste0("sample", 1:6)

# Replace some gene names with real human gene symbols for WebGestaltR
# (WebGestaltR needs recognizable gene symbols)
real_genes <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN", "RB1", "AKT1",
                "MTOR", "PIK3CA", "KRAS", "BRAF", "NRAS", "MAPK1", "MAPK3",
                "JAK2", "STAT3", "CDK4", "CDK6", "CCND1", "CDKN2A",
                "BCL2", "BAX", "CASP3", "CASP9", "TNF", "IL6", "IL1B",
                "TGFB1", "VEGFA", "FGF2", "NOTCH1", "WNT1", "SHH",
                "GAPDH", "ACTB", "TUBB", "HSP90AA1", "HSPA1A", "SOD1",
                "CAT", "GPX1", "NRF2", "KEAP1", "HIF1A", "EPAS1",
                "NFE2L2", "SOX2", "OCT4", "NANOG", "KLF4")
n_real <- min(length(real_genes), n_genes)
rownames(count_matrix)[seq_len(n_real)] <- real_genes[seq_len(n_real)]

col_data <- data.frame(
  condition = factor(rep(c("control", "treated"), each = 3)),
  row.names = colnames(count_matrix)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ condition
)

# Run minimal DESeq2 (fast for 500 genes)
dds <- DESeq(dds, quiet = TRUE)
resF <- results(dds)

save(resF, file = output_path)
cat("Created fixture:", output_path, "\n")
REOF

echo "Creating mock DESeq2 fixture..."
Rscript "${SMOKE_ROOT}/create_fixture.R" "${FIXTURE_RDA}"

if [[ ! -f "${FIXTURE_RDA}" ]]; then
  echo "FAIL: Could not create fixture .rda file"
  exit 1
fi

# ── Run GSEA in smoke mode ───────────────────────────────────────────────────

echo ""
echo "Running GSEA analysis in smoke mode..."
Rscript "${SCRIPT_DIR}/gsea_analysis.R" \
  --deg-results-dir "${FIXTURE_DIR}" \
  --phenotype "tot_adverse_exp" \
  --sex "Female" \
  --output-dir "${RESULTS_DIR}" \
  --smoke

# ── Verify outputs ───────────────────────────────────────────────────────────

echo ""
echo "── Verifying outputs ──"

PASS=true

check_file() {
  if [[ -f "$1" ]]; then
    echo "  OK: $(basename "$1")"
  else
    echo "  MISSING: $(basename "$1")"
    PASS=false
  fi
}

check_file "${RESULTS_DIR}/gsea_summary.csv"
check_file "${RESULTS_DIR}/Exc_ranked_genes.csv"

# Check that at least one RDS was produced
RDS_COUNT=$(find "${RESULTS_DIR}" -name "Exc_*.rds" 2>/dev/null | wc -l)
if [[ "${RDS_COUNT}" -gt 0 ]]; then
  echo "  OK: ${RDS_COUNT} per-database RDS file(s) found"
else
  echo "  WARN: No per-database RDS files (may be expected if no significant results)"
fi

echo ""
if ${PASS}; then
  echo "ACE GSEA smoke test PASSED"
else
  echo "ACE GSEA smoke test FAILED"
  exit 1
fi
