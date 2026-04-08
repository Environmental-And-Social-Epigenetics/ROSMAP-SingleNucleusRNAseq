#!/bin/bash
# generate_report.sh -- Use Claude Code agent to generate a LaTeX report from DEG results
#
# This script reads the machine-readable summaries (deg_summary.json, go_summary.json),
# lists available figures, and calls the Claude CLI to generate a complete LaTeX report.
#
# Usage:
#   bash generate_report.sh \
#     --results-dir results_batch \
#     --figures-dir results_batch/figures \
#     --go-dir results_batch/go_enrichment \
#     --integration "Derived-batch (flowcell)" \
#     --author "Mahmoud Abdelmoneum" \
#     --output-dir report/

set -euo pipefail

# ── Parse arguments ──────────────────────────────────────────────────────────
RESULTS_DIR=""
FIGURES_DIR=""
GO_DIR=""
INTEGRATION=""
AUTHOR="Mahmoud Abdelmoneum"
OUTPUT_DIR="report"
AFFILIATION="Tsai Laboratory, MIT"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir)   RESULTS_DIR="$2"; shift 2 ;;
    --figures-dir)   FIGURES_DIR="$2"; shift 2 ;;
    --go-dir)        GO_DIR="$2"; shift 2 ;;
    --integration)   INTEGRATION="$2"; shift 2 ;;
    --author)        AUTHOR="$2"; shift 2 ;;
    --affiliation)   AFFILIATION="$2"; shift 2 ;;
    --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ -z "$RESULTS_DIR" ]]; then
  echo "ERROR: --results-dir required"
  exit 1
fi

if [[ -z "$FIGURES_DIR" ]]; then
  FIGURES_DIR="${RESULTS_DIR}/figures"
fi
if [[ -z "$GO_DIR" ]]; then
  GO_DIR="${RESULTS_DIR}/go_enrichment"
fi
if [[ -z "$INTEGRATION" ]]; then
  INTEGRATION="$(basename "$RESULTS_DIR" | sed 's/^results_//')"
fi

SUMMARY_JSON="${RESULTS_DIR}/deg_summary.json"
GO_JSON="${GO_DIR}/go_summary.json"

mkdir -p "$OUTPUT_DIR"

# ── Validate inputs ──────────────────────────────────────────────────────────
for f in "$SUMMARY_JSON"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Required file not found: $f"
    echo "Run summarize_results.R first."
    exit 1
  fi
done

FIGURE_LIST=$(ls "$FIGURES_DIR"/*.pdf 2>/dev/null | xargs -I{} basename {} || echo "none")
GO_DATA="No GO enrichment data available."
if [[ -f "$GO_JSON" ]]; then
  GO_DATA=$(cat "$GO_JSON")
fi

DATE=$(date +"%B %Y")
INT_SAFE=$(echo "$INTEGRATION" | tr ' /()' '____' | tr '[:upper:]' '[:lower:]')
TEX_FILE="${OUTPUT_DIR}/${INT_SAFE}_report_draft.tex"
PDF_FILE="${OUTPUT_DIR}/${INT_SAFE}_report_draft.pdf"

echo "=== Report Generation ==="
echo "Integration: $INTEGRATION"
echo "Author:      $AUTHOR"
echo "Results:     $RESULTS_DIR"
echo "Figures:     $FIGURES_DIR"
echo "GO data:     $GO_DIR"
echo "Output:      $TEX_FILE"
echo ""

# ── Locate Claude CLI ────────────────────────────────────────────────────────
CLAUDE_BIN=$(which claude 2>/dev/null || echo "")
if [[ -z "$CLAUDE_BIN" ]]; then
  # Try known locations
  for candidate in \
    /orcd/home/002/mabdel03/conda_envs/consortium/bin/claude \
    "$HOME/.local/bin/claude" \
    "$HOME/.claude/bin/claude"; do
    if [[ -x "$candidate" ]]; then
      CLAUDE_BIN="$candidate"
      break
    fi
  done
fi

if [[ -z "$CLAUDE_BIN" ]]; then
  echo "ERROR: Claude CLI not found. Install with: npm install -g @anthropic-ai/claude-code"
  exit 1
fi

echo "Using Claude CLI: $CLAUDE_BIN"

# ── Build the prompt ──────────────────────────────────────────────────────────
# Truncate JSON to avoid exceeding context limits (keep first 8000 chars)
DEG_DATA=$(head -c 8000 "$SUMMARY_JSON")

PROMPT=$(cat <<'PROMPT_END'
You are a scientific manuscript writer specializing in computational neuroscience and single-cell transcriptomics.

Generate a COMPLETE, compilable LaTeX document body (everything between \begin{document} and \end{document}, excluding those tags) for a DEG analysis report. The report should be publication-quality with detailed biological interpretation.

FORMATTING REQUIREMENTS:
- Use parskip (no paragraph indentation, line-separated paragraphs)
- Do NOT use em dashes anywhere. Use commas, semicolons, or rewrite sentences instead.
- Use \textbf{} for emphasis, not \emph{}
- Reference figures with \ref{} and include them with \includegraphics
- Use \graphicspath for figure paths relative to the figures directory
- Include proper \section, \subsection structure
- Use booktabs-style tables (\toprule, \midrule, \bottomrule)
- Float figures with [H] placement

REQUIRED SECTIONS:
1. Introduction - ACE background, ROSMAP cohort, study motivation
2. Methods - snRNA-seq processing, pseudobulk DESeq2 model (~ age_death + pmi + phenotype + niareagansc), integration method, cell types analyzed
3. Results Overview - Total DEGs by phenotype x sex, sex dimorphism discussion
4. Cell-Type-Specific Results - Table of DEGs per cell type, heatmaps, volcano plots
5. Functional Enrichment - GO analysis interpretation for top cell types
6. PV+ Interneuron Analysis - Dedicated section on PV+ basket cell vulnerability (if data available)
7. Discussion - Biological interpretation, sex dimorphism, microglial/neuronal implications
8. Conclusion
9. Methods Reproducibility

PROMPT_END
)

# Append the actual data
PROMPT="${PROMPT}

ANALYSIS METADATA:
- Integration method: ${INTEGRATION}
- Author: ${AUTHOR}
- Affiliation: ${AFFILIATION}
- Date: ${DATE}

DEG SUMMARY DATA (JSON):
${DEG_DATA}

GO ENRICHMENT DATA (JSON):
${GO_DATA}

AVAILABLE FIGURES (use these exact filenames with \\includegraphics):
${FIGURE_LIST}

The graphicspath should be set to: {${FIGURES_DIR}/}

Output ONLY the LaTeX code, no markdown fences or explanation."

# ── Call Claude ───────────────────────────────────────────────────────────────
echo "Calling Claude agent to generate LaTeX..."

BODY=$("$CLAUDE_BIN" --print --model claude-sonnet-4-6 --max-turns 1 -p "$PROMPT" 2>/dev/null)

if [[ -z "$BODY" ]]; then
  echo "ERROR: Claude returned empty output"
  exit 1
fi

# ── Assemble the full LaTeX document ──────────────────────────────────────────
# Read template and substitute
TEMPLATE_FILE="$(dirname "$0")/report/template.tex"
if [[ ! -f "$TEMPLATE_FILE" ]]; then
  TEMPLATE_FILE="$(cd "$(dirname "$0")" && pwd)/report/template.tex"
fi

if [[ -f "$TEMPLATE_FILE" ]]; then
  # Use the template: replace placeholders
  sed \
    -e "s|%%% __GRAPHICSPATH__ %%%|\\\\graphicspath{{${FIGURES_DIR}/}}|" \
    -e "s|%%% __TITLE__ %%%|\\\\title{\\\\textbf{Cell-Type-Specific Transcriptomic Effects of Adverse Childhood Experiences in the Human Prefrontal Cortex} \\\\\\\\\\\\[0.5em] \\\\large ${INTEGRATION} Integration}|" \
    -e "s|%%% __AUTHOR__ %%%|\\\\author{${AUTHOR} \\\\\\\\ ${AFFILIATION}}|" \
    -e "s|%%% __DATE__ %%%|\\\\date{${DATE}}|" \
    "$TEMPLATE_FILE" > "$TEX_FILE.tmp"

  # Insert body before \end{document}
  sed -i "/%%% __BODY__ %%%/r /dev/stdin" "$TEX_FILE.tmp" <<< "$BODY"
  sed -i '/%%% __BODY__ %%%/d' "$TEX_FILE.tmp"
  mv "$TEX_FILE.tmp" "$TEX_FILE"
else
  # No template: write the agent's output directly (it should be a complete document)
  echo "$BODY" > "$TEX_FILE"
fi

echo "LaTeX written to: $TEX_FILE"

# ── Compile PDF ───────────────────────────────────────────────────────────────
XELATEX=$(which xelatex 2>/dev/null || echo "$HOME/.local/bin/xelatex")
if [[ -x "$XELATEX" ]]; then
  echo "Compiling PDF..."
  cd "$OUTPUT_DIR"
  TEX_BASE=$(basename "$TEX_FILE")
  "$XELATEX" -interaction=nonstopmode "$TEX_BASE" > /dev/null 2>&1 || true
  "$XELATEX" -interaction=nonstopmode "$TEX_BASE" > /dev/null 2>&1 || true
  cd - > /dev/null

  if [[ -f "$PDF_FILE" ]]; then
    echo "PDF compiled: $PDF_FILE"
  else
    echo "WARNING: PDF compilation failed. Check $TEX_FILE for errors."
  fi
else
  echo "WARNING: xelatex not found. LaTeX written but not compiled."
fi

echo ""
echo "Report generation complete."
echo "  LaTeX: $TEX_FILE"
echo "  PDF:   $PDF_FILE (if compiled)"
echo ""
echo "NOTE: This is a DRAFT. Review and edit before treating as final."
