#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}"

bash "${SCRIPT_DIR}/run_functional_concordance.sh" "$@"
