#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellbender_config.sh"

total=0
if [[ -f "${BATCH_ASSIGNMENTS_CSV}" ]]; then
    total=$(($(wc -l < "${BATCH_ASSIGNMENTS_CSV}") - 1))
fi

completed=0
failed=0
if [[ -f "${TRACKING_DIR}/cellbender_completed.txt" ]]; then
    completed=$(wc -l < "${TRACKING_DIR}/cellbender_completed.txt")
fi
if [[ -f "${TRACKING_DIR}/cellbender_failed.txt" ]]; then
    failed=$(wc -l < "${TRACKING_DIR}/cellbender_failed.txt")
fi

remaining=$(( total - completed - failed ))
if [[ ${remaining} -lt 0 ]]; then
    remaining=0
fi

echo "CellBender status:"
echo "  Total samples:     ${total}"
echo "  Completed:         ${completed}"
echo "  Failed:            ${failed}"
echo "  Remaining:         ${remaining}"
