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

for tool in "$PYTHON" meme fimo; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Required command not found: $tool" >&2
        exit 1
    fi
done

"$PYTHON" -c 'import pandas, numpy, Bio, openpyxl'

for input in \
    results/operon_prediction_ll/likehood.json \
    data/S_roseum/CDS_with_protein.csv \
    data/S_roseum/genomic.fna \
    data/S_lividans/genome.fasta \
    data/S_lividans/TSS_table.xlsx; do
    if [[ ! -s "$input" ]]; then
        echo "Required input missing or empty: $input" >&2
        exit 1
    fi
done

mkdir -p results/promoter_prediction

echo "[1/6] Predicting operons"
"$PYTHON" scripts/operon_prediction/ll_using_analysis.py

echo "[2/6] Extracting promoter training and upstream sequences"
"$PYTHON" scripts/promoter_prediction/extract_upstream.py

echo "[3/6] Discovering motifs with MEME"
bash scripts/promoter_prediction/run_meme.sh

echo "[4/6] Scanning upstream sequences with FIMO"
bash scripts/promoter_prediction/run_fimo.sh

echo "[5/6] Pairing promoter elements"
"$PYTHON" scripts/promoter_prediction/process_fimo.py

echo "[6/6] Exporting operons with promoters"
"$PYTHON" scripts/operon_prediction/operon_prediction.py

echo "Analysis complete: $PROJECT_ROOT/results/operons_with_promoters.csv"
