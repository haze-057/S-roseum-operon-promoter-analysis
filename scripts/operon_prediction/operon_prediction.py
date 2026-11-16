
import csv
import json
from collections import defaultdict
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[2]


def operons_with_promoters(operon_path: Path, promoter_path: Path,
                          replicon: str, strand: str) -> list:

    with operon_path.open(encoding="utf-8") as f:
        operons = json.load(f)

    with promoter_path.open(encoding="utf-8") as f:
        promoters = json.load(f)

    promoter_groups = defaultdict(list)

    for hit35, hit10 in promoters[1:]:
        operon_id = str(hit35[2])

        if str(hit10[2]) != operon_id:
            raise ValueError("Promoter elements belong to different operons")

        promoter_groups[operon_id].append((hit35, hit10))

    results = []

    for i, operon in enumerate(operons):
        operon_id = str(i)
        genes = [gene[4] for gene in operon]

        operon_start = min(gene[1] for gene in operon)
        operon_end = max(gene[2] for gene in operon)

        matches = promoter_groups.get(operon_id, [])

        if not matches:
            matches = [(None, None)]

        for hit35, hit10 in matches:

            if hit35:
                if strand == "plus":
                    genomic35_start = operon_start - 300 + hit35[3]
                    genomic35_end = operon_start - 300 + hit35[4]
                    genomic10_start = operon_start - 300 + hit10[3]
                    genomic10_end = operon_start - 300 + hit10[4]
                else:
                    genomic35_start = operon_end + 301 - hit35[4]
                    genomic35_end = operon_end + 301 - hit35[3]
                    genomic10_start = operon_end + 301 - hit10[4]
                    genomic10_end = operon_end + 301 - hit10[3]
            else:
                genomic35_start = genomic35_end = None
                genomic10_start = genomic10_end = None

            results.append({
                "replicon": replicon,
                "strand": strand,
                "operon_id": operon_id,
                "genes": ", ".join(map(str, genes)),
                "locus_tags": ", ".join(str(gene[5]) for gene in operon),
                "gene_count": len(genes),
                "operon_start": operon_start,
                "operon_end": operon_end,
                "minus35_start": hit35[3] if hit35 else None,
                "minus10_start": hit10[3] if hit10 else None,
                "minus35_genomic_start": genomic35_start,
                "minus35_genomic_end": genomic35_end,
                "minus10_genomic_start": genomic10_start,
                "minus10_genomic_end": genomic10_end,
                "distance": hit10[3] - hit35[3] if hit35 else None,
                "pvalue35": hit35[7] if hit35 else None,
                "pvalue10": hit10[7] if hit10 else None,
            })

    return results


results = []

for replicon in ("chr", "pl"):
    for strand in ("plus", "minus"):
        results.extend(
            operons_with_promoters(
                BASE_DIR / "results/operon_prediction_ll"
                / f"{replicon}_{strand}_distances_ll.json",

                BASE_DIR / "results/promoter_prediction"
                / f"pr_{replicon}_{strand}",

                replicon,
                strand
            )
        )

output_path = BASE_DIR / "results/operons_with_promoters.csv"

with output_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print(f"Saved {len(results)} rows to {output_path}")
