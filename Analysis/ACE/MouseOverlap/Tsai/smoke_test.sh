#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

set +u
activate_env "${NEBULA_ENV}"
set -u

TMP_ROOT="${TMPDIR:-/tmp}/ace_mouse_overlap_smoke_$$"
MOUSE_DIR="${TMP_ROOT}/mouse"
OUT_DIR="${TMP_ROOT}/out"
mkdir -p "${MOUSE_DIR}" "${OUT_DIR}/figures"

cleanup() {
  rm -rf "${TMP_ROOT}"
}
trap cleanup EXIT

cat > "${MOUSE_DIR}/PairedEnd_GeneLevel_PFC_NeuN_NGHvLNB_Dark_DESEQ2output.csv" <<'CSV'
,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,gene.name
ENSMUSG00000000001,100,-2.0,0.2,-10,1e-8,0.001,M1
ENSMUSG00000000002,100,1.0,0.2,5,1e-5,0.01,M2
ENSMUSG00000000003,100,-1.0,0.2,-5,1e-5,0.2,M3
ENSMUSG00000000004,100,-1.0,0.2,-5,1e-5,0.01,M4
CSV

cat > "${TMP_ROOT}/HOM_MouseHumanSequence.rpt" <<'TSV'
DB Class Key	Common Organism Name	NCBI Taxon ID	Symbol	EntrezGene ID	Mouse MGI ID	HGNC ID	OMIM Gene ID	Genetic Location	Genome Coordinates (mouse: GRCm39 human: GRCh38)	Nucleotide RefSeq IDs	Protein RefSeq IDs	SWISS_PROT IDs
1	mouse, laboratory	10090	M1	1	MGI:1							
1	human	9606	H1	11		HGNC:1						
2	mouse, laboratory	10090	M2	2	MGI:2							
2	human	9606	H2	12		HGNC:2						
3	mouse, laboratory	10090	M3	3	MGI:3							
3	human	9606	H3	13		HGNC:3						
4	mouse, laboratory	10090	M4	4	MGI:4							
4	human	9606	H4	14		HGNC:4						
5	mouse, laboratory	10090	M5A	5	MGI:5							
5	mouse, laboratory	10090	M5B	6	MGI:6							
5	human	9606	H5	15		HGNC:5						
TSV

cat > "${TMP_ROOT}/human_deg_flat.csv" <<'CSV'
integration,phenotype,cell_type,sex,gene_symbol,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,human_sig,human_direction
derived_batch,tot_adverse_exp,Exc,Male,H1,100,1.0,0.1,10,1e-8,0.001,TRUE,ACE_up
derived_batch,tot_adverse_exp,Exc,Male,H2,100,-1.0,0.1,-10,1e-8,0.001,TRUE,ACE_down
derived_batch,tot_adverse_exp,Exc,Male,H3,100,0.5,0.1,5,1e-5,0.01,TRUE,ACE_up
derived_batch,tot_adverse_exp,Exc,Male,H4,100,-0.5,0.1,-5,1e-5,0.2,FALSE,ACE_down
CSV

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/normalize_mouse_deg.py" \
  --mouse-dir "${MOUSE_DIR}" \
  --output "${OUT_DIR}/mouse_deg_normalized.csv"

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/prepare_orthologs.py" \
  --source "${TMP_ROOT}/HOM_MouseHumanSequence.rpt" \
  --raw-output "${OUT_DIR}/HOM_MouseHumanSequence.rpt" \
  --output "${OUT_DIR}/ortholog_mapping_used.csv" \
  --metadata-output "${OUT_DIR}/ortholog_mapping_metadata.json"

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/overlap_analysis.py" \
  --mouse-deg "${OUT_DIR}/mouse_deg_normalized.csv" \
  --human-deg "${TMP_ROOT}/human_deg_flat.csv" \
  --ortholog-table "${OUT_DIR}/ortholog_mapping_used.csv" \
  --output-dir "${OUT_DIR}" \
  --figures-dir "${OUT_DIR}/figures"

"${NEBULA_ENV}/bin/python" - <<PY
import pandas as pd
from pathlib import Path
out = Path("${OUT_DIR}")
mouse = pd.read_csv(out / "mouse_deg_normalized.csv")
assert mouse.loc[mouse["mouse_symbol"] == "M1", "mouse_lnb_log2fc"].iloc[0] == 2.0
assert mouse.loc[mouse["mouse_symbol"] == "M2", "mouse_direction"].iloc[0] == "LNB_down"
orth = pd.read_csv(out / "ortholog_mapping_used.csv")
assert int(orth["is_one_to_one"].sum()) == 4
summary = pd.read_csv(out / "overlap_summary.csv")
row = summary[
    (summary["mouse_dataset"] == "NeuN_Dark")
    & (summary["human_cell_type"] == "Exc")
    & (summary["human_sex"] == "Male")
    & (summary["mouse_direction"] == "LNB_up")
    & (summary["human_direction"] == "ACE_up")
].iloc[0]
assert int(row["universe_n"]) == 4
assert int(row["mouse_set_n"]) == 2
assert int(row["human_set_n"]) == 2
assert int(row["overlap_n"]) == 1
genes = pd.read_csv(out / "overlap_genes.csv")
assert "H1" in set(genes["human_symbol"])
print("Smoke assertions passed.")
PY

echo "Smoke test passed."

