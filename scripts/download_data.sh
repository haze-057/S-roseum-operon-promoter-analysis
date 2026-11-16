#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"
DATA_DIR="$PROJECT_ROOT/data"

if [[ -z "${PYTHON:-}" ]]; then
    if [[ -x "$PROJECT_ROOT/.venv/bin/python" ]]; then
        PYTHON="$PROJECT_ROOT/.venv/bin/python"
    else
        PYTHON=python3
    fi
fi
for tool in curl gunzip "$PYTHON"; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Required command not found: $tool" >&2; exit 1; }
done

mkdir -p "$DATA_DIR/S_coelicolor"
mkdir -p "$DATA_DIR/S_lividans"
mkdir -p "$DATA_DIR/S_roseum"

download() {
    local url="$1"
    local output="$2"

    echo "Downloading $(basename "$output")"

    curl -fL --retry 3 \
        --connect-timeout 20 \
        "$url" -o "$output"
}

# Streptomyces coelicolor

download \
    "https://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-11-578/MediaObjects/12864_2010_9973_MOESM2_ESM.CSV" \
    "$DATA_DIR/S_coelicolor/gene_pairs.CSV"

download \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AL645882.2&rettype=gbwithparts&retmode=text" \
    "$DATA_DIR/S_coelicolor/sequence.gb"


# Streptomyces lividans

download \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP009124.1&rettype=fasta&retmode=text" \
    "$DATA_DIR/S_lividans/genome.fasta"

download \
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/6742748/bin/Table_4.XLSX" \
    "$DATA_DIR/S_lividans/TSS_table.xlsx"


# Streptosporangium roseum

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/024/865/GCF_000024865.1_ASM2486v1"
PREFIX="GCF_000024865.1_ASM2486v1"

for filename in genomic.fna cds_from_genomic.fna feature_table.txt; do
    download \
        "$BASE_URL/${PREFIX}_${filename}.gz" \
        "$DATA_DIR/S_roseum/${filename}.gz"

    gunzip -f "$DATA_DIR/S_roseum/${filename}.gz"
done

"$PYTHON" "$PROJECT_ROOT/scripts/prepare_cds.py" \
    "$DATA_DIR/S_roseum/feature_table.txt" "$DATA_DIR/S_roseum/CDS_with_protein.csv"

echo "All source datasets downloaded and CDS table prepared."
