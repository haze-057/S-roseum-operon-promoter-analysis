
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_DIR="$PROJECT_ROOT/data/S_lividans"
OUTPUT_DIR="$PROJECT_ROOT/results/promoter_prediction"

meme "$INPUT_DIR/promoters_minus10.fasta" \
    -dna -mod zoops \
    -minw 6 -maxw 6 \
    -oc "$OUTPUT_DIR/meme_minus10"

meme "$INPUT_DIR/promoters_minus35.fasta" \
    -dna -mod oops \
    -minw 4 -maxw 4 \
    -oc "$OUTPUT_DIR/meme_minus35"
