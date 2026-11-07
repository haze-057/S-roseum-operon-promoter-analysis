import json

import pandas as pd


def distances_search(df: pd.DataFrame, distance: int, strain: str) -> list:
    df = df.astype(object) 
    op = []
    for i in range(df.shape[0]-1):
        if abs(df.iloc[i, df.columns.get_loc('end')] - df.iloc[i+1, df.columns.get_loc('start')]) <= distance and df.iloc[i, df.columns.get_loc('strand')] == strain and df.iloc[i+1, df.columns.get_loc('strand')] == strain:
            op.append(tuple(df.iloc[i+1]))
        else:
            if len(op) > 1:
                yield op
            op = []
            op.append(tuple(df.iloc[i]))

cds_with_protein = pd.read_csv("CDS_with_protein.csv")

cds_sel = cds_with_protein[['seq_type', 'start', 'end', 'strand', 'name']].copy()

chromosome_proteins = cds_sel[cds_sel['seq_type'] == 'chromosome']
plasmid_proteins = cds_sel[cds_sel['seq_type'] == 'plasmid']

distance=50

with open("jsons/chr_plus_distances.json", "w") as f:
    json.dump(list(distances_search(df=chromosome_proteins, distance=distance, strain='+')), f)

with open("jsons/chr_minus_distances.json", "w") as f:
    json.dump(list(distances_search(df=chromosome_proteins, distance=distance, strain='-')), f)

with open("jsons/pl_plus_distances.json", "w") as f:
    json.dump(list(distances_search(df=plasmid_proteins, distance=distance, strain='+')), f)

with open("jsons/pl_minus_distances.json", "w") as f:
    json.dump(list(distances_search(df=plasmid_proteins, distance=distance, strain='-')), f)

