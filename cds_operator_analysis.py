import pandas as pd


def distances_search(df: pd.DataFrame, distance: int) -> list:
    op = []
    for i in range(df.shape[0]-1):
        if abs(df.iloc[i, df.columns.get_loc('end')] - df.iloc[i+1, df.columns.get_loc('start')]) <= distance and df.iloc[i, df.columns.get_loc('strand')] == df.iloc[i+1, df.columns.get_loc('strand')]:
            op.append(tuple(df.iloc[i+1]))
        else:
            if len(op) != 1:
                yield op
            op = []
            op.append(tuple(df.iloc[i]))

cds_with_protein = pd.read_csv("CDS_with_protein.csv")

cds_sel = cds_with_protein[['seq_type', 'start', 'end', 'strand', 'name']].copy()

chromosome_proteins = cds_sel[cds_sel['seq_type'] == 'chromosome']
plasmid_proteins = cds_sel[cds_sel['seq_type'] == 'plasmid']

distance=50

chromosome_plus = list(distances_search(df=chromosome_proteins, distance=distance))
