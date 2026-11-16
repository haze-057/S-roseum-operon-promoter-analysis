import json
from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[2]


def promoters_on_seq(el_10_file_name: str, el_35_file_name: str, output_name: str):
    el_10_table = pd.read_csv(
        BASE_DIR / el_10_file_name, sep="\t", comment="#", dtype={"sequence_name": str}
    )
    el_35_table = pd.read_csv(
        BASE_DIR / el_35_file_name, sep="\t", comment="#", dtype={"sequence_name": str}
    )
    
    m = el_10_table.shape[0]

    el_10_table = el_10_table.sort_values(
        ['sequence_name', 'start']
    ).reset_index(drop=True)

    el_35_table = el_35_table.sort_values(
        ['sequence_name', 'start']
    ).reset_index(drop=True)
    
    promoters_el_10_35 = [["motif_id",
      "motif_alt_id",
      "sequence_name",
      "start",
      "stop",
      "strand",
      "score",
      "p-value",
      "q-value", 
      "matched_sequence"]]

    for seq_name, hits35 in el_35_table.groupby("sequence_name"):
        hits10 = el_10_table[
            el_10_table["sequence_name"] == seq_name
        ].reset_index(drop=True)

        m = hits10.shape[0]
        lo = hi = 0

        for i in range(hits35.shape[0]):
            while lo < m and hits10.iloc[lo]['start'] < hits35.iloc[i]['start'] + 20:
                lo += 1
            while hi < m and hits10.iloc[hi]['start'] <= hits35.iloc[i]['start'] + 26:
                hi += 1
            for j in range(lo, hi):
                promoters_el_10_35.append([
                    list(hits35.iloc[i].to_dict().values()),
                    list(hits10.iloc[j].to_dict().values()),
                ])

    print(len(promoters_el_10_35)-1)

    with (BASE_DIR / output_name).open("w", encoding="utf-8") as output:
        json.dump(promoters_el_10_35, output, indent=2)

for a in ['chr', 'pl']:
    for b in ['plus', 'minus']:
        promoters_on_seq(f"results/promoter_prediction/fimo_{a}_{b}_minus10/fimo.tsv", f"results/promoter_prediction/fimo_{a}_{b}_minus35/fimo.tsv", f"results/promoter_prediction/pr_{a}_{b}")
