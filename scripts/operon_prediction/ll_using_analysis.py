import json
from pathlib import Path

import numpy as np
import pandas as pd


def get_ll(ll_list: list, dist: int, w_bins: int = 10) -> float:
    if dist >= ll_list[0] and dist <= ll_list[1]:
        return ll_list[-1][(dist - ll_list[0]) // w_bins]
    return float("-inf")

def distances_search(pr_list: pd.DataFrame, ll_list: list, strand: str, ll: float = -0.03) -> list:
    op = []
    for i in range(pr_list.shape[0] - 1):
        if pr_list.iloc[i]['strand'] == strand and pr_list.iloc[i + 1]['strand'] == strand:
            dist = pr_list.iloc[i + 1]['start'] - pr_list.iloc[i]['end'] - 1
            if get_ll(ll_list, dist) > ll:
                if not op:
                    op.append(tuple(pr_list.iloc[i]))
                op.append(tuple(pr_list.iloc[i+1]))
            else:
                if len(op) > 1:
                    yield op
                op = []
    if len(op) > 1:
        yield op

cds_with_protein = pd.read_csv("data/S_roseum/CDS_with_protein.csv")

cds_sel = cds_with_protein[['seq_type', 'start', 'end', 'strand', 'name', 'locus_tag']].copy()

chromosome_proteins = cds_sel[cds_sel['seq_type'] == 'chromosome']
plasmid_proteins = cds_sel[cds_sel['seq_type'] == 'plasmid']

ll_path = Path(__file__).resolve().parents[2] / "results" / "operon_prediction_ll" 

with open(ll_path / "likehood.json", "r") as f:
    ll_results = json.load(f)

with open(ll_path / "chr_plus_distances_ll.json", "w") as f:
    json.dump(list(distances_search(pr_list=chromosome_proteins, ll_list=ll_results, strand='+')), f, default=lambda value: value.item())

with open(ll_path / "chr_minus_distances_ll.json", "w") as f:
    json.dump(list(distances_search(pr_list=chromosome_proteins, ll_list=ll_results, strand='-')), f, default=lambda value: value.item())

with open(ll_path / "pl_plus_distances_ll.json", "w") as f:
    json.dump(list(distances_search(pr_list=plasmid_proteins, ll_list=ll_results, strand='+')), f, default=lambda value: value.item())

with open(ll_path / "pl_minus_distances_ll.json", "w") as f:
    json.dump(list(distances_search(pr_list=plasmid_proteins, ll_list=ll_results, strand='-')), f, default=lambda value: value.item())

groups = list(distances_search(
    pr_list=chromosome_proteins,
    ll_list=ll_results,
    strand='+',
))

print(np.mean([len(group) for group in groups]) if groups else float("nan"))
