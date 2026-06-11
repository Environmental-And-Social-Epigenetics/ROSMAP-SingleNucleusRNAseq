#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

SMOKE_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/Smoke/Tsai/MouseFunctionalTotAdverseMale"
INPUT_DIR="${SMOKE_ROOT}/inputs"
HUMAN_DIR="${INPUT_DIR}/human_gsea/Male"
OUTPUT_DIR="${SMOKE_ROOT}/outputs"
FIGURES_DIR="${OUTPUT_DIR}/figures"

rm -rf "${SMOKE_ROOT}"
mkdir -p "${INPUT_DIR}" "${HUMAN_DIR}" "${FIGURES_DIR}"

cat > "${INPUT_DIR}/mouse.csv" <<'CSV'
mouse_dataset,mouse_population,condition,primary_human_target,source_file,mouse_ensembl_id,mouse_symbol,mouse_baseMean,mouse_original_log2FoldChange,mouse_lnb_log2fc,mouse_lfcSE,mouse_stat,mouse_pvalue,mouse_padj,mouse_sig,mouse_direction
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0001,M1,100,-2,2,0.1,-6,1e-8,0.001,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0002,M2,100,-2,2,0.1,-5,1e-7,0.002,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0003,M3,100,-2,2,0.1,-4,1e-6,0.003,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0004,M4,100,-2,2,0.1,-3,1e-5,0.004,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0005,M5,100,-2,2,0.1,-2.8,1e-5,0.005,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0006,M6,100,-2,2,0.1,-2.6,1e-5,0.006,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0007,M7,100,-2,2,0.1,-2.4,1e-5,0.007,TRUE,LNB_up
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0008,M8,100,1,-1,0.1,1,0.1,0.5,FALSE,LNB_down
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0009,M9,100,1,-1,0.1,1.2,0.1,0.5,FALSE,LNB_down
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0010,M10,100,1,-1,0.1,1.4,0.1,0.5,FALSE,LNB_down
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0011,M11,100,1,-1,0.1,1.6,0.1,0.5,FALSE,LNB_down
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0012,M12,100,1,-1,0.1,1.8,0.1,0.5,FALSE,LNB_down
NeuN_Dark,NeuN,Dark,Exc,synthetic,ENSMUSG0013,M13,100,-3,3,0.1,-8,1e-9,0.001,TRUE,LNB_up
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0001,M1,100,-1,1,0.1,-3,1e-4,0.02,TRUE,LNB_up
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0002,M2,100,-1,1,0.1,-2.5,1e-3,0.03,TRUE,LNB_up
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0003,M3,100,-1,1,0.1,-2,1e-2,0.04,TRUE,LNB_up
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0004,M4,100,-1,1,0.1,-1.8,1e-2,0.05,FALSE,LNB_up
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0005,M5,100,1,-1,0.1,1.2,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0006,M6,100,1,-1,0.1,1.4,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0007,M7,100,1,-1,0.1,1.6,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0008,M8,100,1,-1,0.1,1.8,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0009,M9,100,1,-1,0.1,2,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0010,M10,100,1,-1,0.1,2.2,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0011,M11,100,1,-1,0.1,2.4,0.1,0.5,FALSE,LNB_down
NeuN_Light,NeuN,Light,Exc,synthetic,ENSMUSG0012,M12,100,1,-1,0.1,2.6,0.1,0.5,FALSE,LNB_down
CSV

cat > "${INPUT_DIR}/ortholog.csv" <<'CSV'
ortholog_source,db_class_key,mouse_symbol,human_symbol,mouse_entrez,human_entrez,mouse_mgi_id,hgnc_id,n_mouse_in_class,n_human_in_class,is_one_to_one
synthetic,1,M1,H1,1,101,MGI:1,HGNC:1,1,1,TRUE
synthetic,2,M2,H2,2,102,MGI:2,HGNC:2,1,1,TRUE
synthetic,3,M3,H3,3,103,MGI:3,HGNC:3,1,1,TRUE
synthetic,4,M4,H4,4,104,MGI:4,HGNC:4,1,1,TRUE
synthetic,5,M5,H5,5,105,MGI:5,HGNC:5,1,1,TRUE
synthetic,6,M6,H6,6,106,MGI:6,HGNC:6,1,1,TRUE
synthetic,7,M7,H7,7,107,MGI:7,HGNC:7,1,1,TRUE
synthetic,8,M8,H8,8,108,MGI:8,HGNC:8,1,1,TRUE
synthetic,9,M9,H9,9,109,MGI:9,HGNC:9,1,1,TRUE
synthetic,10,M10,H10,10,110,MGI:10,HGNC:10,1,1,TRUE
synthetic,11,M11,H11,11,111,MGI:11,HGNC:11,1,1,TRUE
synthetic,12,M12,H12,12,112,MGI:12,HGNC:12,1,1,TRUE
synthetic,13,M13,H13,13,113,MGI:13,HGNC:13,2,1,FALSE
CSV

cat > "${INPUT_DIR}/toy.gmt" <<'GMT'
TOY_CONCORDANT_UP	Toy concordant up	H1	H2	H3	H4
TOY_DISCORDANT	Toy discordant	H5	H6	H7
TOY_BACKGROUND	Toy background	H8	H9	H10	H11	H12
GMT

set +u
activate_env "${GSEA_ANALYSIS_ENV}"
set -u

Rscript - <<RS
dir.create("${HUMAN_DIR}", recursive = TRUE, showWarnings = FALSE)
mk <- function(nes_up = 2, nes_dis = -2) {
  data.frame(
    geneSet = c("TOY_CONCORDANT_UP", "TOY_DISCORDANT", "TOY_BACKGROUND"),
    description = c("Toy concordant up", "Toy discordant", "Toy background"),
    normalizedEnrichmentScore = c(nes_up, nes_dis, -0.2),
    pValue = c(0.001, 0.002, 0.5),
    FDR = c(0.01, 0.02, 0.8),
    size = c(4, 3, 5),
    leadingEdgeNum = c(4, 3, 2),
    userId = c("H1;H2;H3;H4", "H5;H6;H7", "H8;H9"),
    stringsAsFactors = FALSE
  )
}
saveRDS(mk(2, -2), file.path("${HUMAN_DIR}", "Ex-L2_3_custom.rds"))
saveRDS(mk(1.5, -1.5), file.path("${HUMAN_DIR}", "Ex-L4_custom.rds"))
saveRDS(mk(-2, 2), file.path("${HUMAN_DIR}", "Inh_custom.rds"))
RS

Rscript "${SCRIPT_DIR}/functional_tot_adverse_male.R" \
  --mouse-deg "${INPUT_DIR}/mouse.csv" \
  --ortholog-table "${INPUT_DIR}/ortholog.csv" \
  --human-gsea-dir "${HUMAN_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --figures-dir "${FIGURES_DIR}" \
  --custom-gmt "${INPUT_DIR}/toy.gmt" \
  --min-size 3 \
  --max-size 20 \
  --permutations 100 \
  --fdr-threshold 1

python - <<PY
from pathlib import Path
import pandas as pd

out = Path("${OUTPUT_DIR}")
ranked = pd.read_csv(out / "mouse_ranked_genes.csv")
mouse_gsea = pd.read_csv(out / "mouse_gsea_results.csv")
ex = pd.read_csv(out / "human_ex_aggregate_reference.csv")
summary = pd.read_csv(out / "functional_concordance_summary.csv")
paths = pd.read_csv(out / "functional_concordance_pathways.csv")

assert set(["NeuN_Dark", "NeuN_Light"]).issubset(set(ranked["mouse_dataset"])), "dark/light datasets were not kept separate"
assert "H13" not in set(ranked["human_symbol"]), "non-one-to-one ortholog was not filtered"
assert ranked.loc[(ranked["mouse_dataset"] == "NeuN_Dark") & (ranked["human_symbol"] == "H1"), "rank"].iloc[0] > 0, "mouse LNB rank sign flip failed"
assert len(ranked[ranked["mouse_dataset"] == "NeuN_Dark"]) == 12, "not all finite one-to-one ranked genes were retained"

toy_up = ex[ex["geneSet"] == "TOY_CONCORDANT_UP"].iloc[0]
assert toy_up["subtype_support_n"] == 2, "Ex aggregate support count is wrong"
assert toy_up["NES"] > 0, "Ex aggregate median NES should be positive"
assert toy_up["subtype_sign_consistency"] == 1.0, "Ex aggregate sign consistency is wrong"

assert (summary["comparison_role"] == "primary").any(), "primary comparison missing"
assert "custom" in set(mouse_gsea["library"]), "custom smoke GSEA library missing"

p_up = paths[(paths["comparison_role"] == "primary") & (paths["geneSet"] == "TOY_CONCORDANT_UP")].iloc[0]
assert bool(p_up["is_concordant"]), "toy concordant pathway not classified as concordant"
p_dis = paths[(paths["comparison_role"] == "primary") & (paths["geneSet"] == "TOY_DISCORDANT")].iloc[0]
assert bool(p_dis["is_discordant"]), "toy discordant pathway not classified as discordant"

print("Functional total-adverse male smoke test PASSED")
PY
