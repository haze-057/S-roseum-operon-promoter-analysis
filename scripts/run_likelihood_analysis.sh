#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

if [[ -z "${PYTHON:-}" ]]; then
    if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
        PYTHON="$PROJECT_ROOT/.venv/bin/python"
    else
        PYTHON=python3
    fi
fi

if ! command -v "$PYTHON" >/dev/null 2>&1; then
    echo "Required command not found: $PYTHON" >&2
    exit 1
fi

"$PYTHON" -c 'import pandas, numpy, Bio, matplotlib'

for input in \
    data/S_coelicolor/sequence.gb \
    data/S_coelicolor/gene_pairs.CSV; do
    if [[ ! -s "$input" ]]; then
        echo "Required input missing or empty: $input" >&2
        exit 1
    fi
done

mkdir -p results/operon_prediction_ll

echo "Calculating inter-CDS distance likelihoods"
"$PYTHON" scripts/likehood_analysis/ll_calc.py

echo "Analysis complete: $PROJECT_ROOT/results/operon_prediction_ll/likehood.json"
