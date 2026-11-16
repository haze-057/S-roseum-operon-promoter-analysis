import argparse
import csv
from pathlib import Path


def prepare_cds(source: Path, output: Path) -> None:
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"# feature", "class", "seq_type", "start", "end", "strand", "name", "locus_tag"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Input is missing required NCBI feature table columns")
        columns = [column for column in reader.fieldnames if column != "GeneID"]
        rows = [
            {column: row[column] for column in columns}
            for row in reader
            if row["# feature"] == "CDS" and row["class"] == "with_protein"
        ]
    if not rows:
        raise ValueError("No protein-coding CDS rows found")
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} CDS rows to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    prepare_cds(args.source, args.output)
