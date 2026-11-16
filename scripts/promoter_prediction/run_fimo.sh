#!/usr/bin/env bash
set -euo pipefail
    
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_DIR="$PROJECT_ROOT/data/S_roseum"
RESULTS_DIR="$PROJECT_ROOT/results/promoter_prediction"

for replicon in chr pl; do
    for strand in plus minus; do
        for motif in minus10 minus35; do

            input="$INPUT_DIR/${replicon}_${strand}_upstream.fasta"
            motif_file="$RESULTS_DIR/meme_${motif}/meme.txt"
            output="$RESULTS_DIR/fimo_${replicon}_${strand}_${motif}"

            echo "Scanning $replicon $strand for $motif"

            fimo \
                --norc \
                --thresh 0.01 \
                --oc "$output" \
                "$motif_file" \
                "$input"

        done
    done
done