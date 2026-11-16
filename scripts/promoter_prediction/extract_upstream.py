import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def export_pr_ll(op_list: list, from_fasta: str, output_path: Path, strand: str = '+'):
    seq = []
    if strand == '+':
        for i, operon in enumerate(op_list):
            seq.append([str(i), from_fasta[max(0, operon[0][1]-300):operon[0][1]]])
    else:
        for i, operon in enumerate(op_list):
            operon_end = max(gene[2] for gene in operon)

            seq.append([
                str(i),
                from_fasta[
                    operon_end:operon_end + 300
                ].reverse_complement()
            ])
    r_sequences = [
        SeqRecord(
            Seq(sequence),
            id=tss_id,
            description=""
        )
        for tss_id, sequence in seq
    ]
    with output_path.open("w", encoding="utf-8") as output:
        SeqIO.write(r_sequences, output, "fasta")

json_dir = Path(__file__).resolve().parents[2] / "results" / "operon_prediction_ll"

with (json_dir / "chr_plus_distances_ll.json").open(encoding="utf-8") as source:
    chr_plus_distances_ll = json.load(source)
with (json_dir / "chr_minus_distances_ll.json").open(encoding="utf-8") as source:
    chr_minus_distances_ll = json.load(source)
with (json_dir / "pl_plus_distances_ll.json").open(encoding="utf-8") as source:
    pl_plus_distances_ll = json.load(source)
with (json_dir / "pl_minus_distances_ll.json").open(encoding="utf-8") as source:
    pl_minus_distances_ll = json.load(source)

sequences = []

for record in SeqIO.parse("data/S_lividans/genome.fasta", "fasta"):
    print(record.id)
    print(len(record.seq))
    print(record.description)
    sequences.append(record.seq)

S_lividans_chromosome_seq = sequences[0].upper()

sequences = []

for record in SeqIO.parse("data/S_roseum/genomic.fna", "fasta"):
    print(record.id)
    print(len(record.seq))
    print(record.description)
    sequences.append(record.seq)

S_roseum_chromosome_seq = sequences[0].upper()
S_roseum_plasmid_seq = sequences[1].upper()

tss_table = pd.read_excel("data/S_lividans/TSS_table.xlsx")

tss_table = tss_table[['Position', 'Intensity', 'Strand', 'Associated gene']]

sequences35 = []
sequences10 = []

for i in range(tss_table.shape[0]):
    if tss_table.iloc[i]['Strand'] == '+':
        sequences35.append([str(i), S_lividans_chromosome_seq[tss_table.iloc[i]['Position']-41:tss_table.iloc[i]['Position']-26]])
        sequences10.append([str(i), S_lividans_chromosome_seq[tss_table.iloc[i]['Position']-21:tss_table.iloc[i]['Position']]])
    else:
        sequences35.append([str(i), S_lividans_chromosome_seq[tss_table.iloc[i]['Position']+25:tss_table.iloc[i]['Position']+40].reverse_complement()])
        sequences10.append([str(i), S_lividans_chromosome_seq[tss_table.iloc[i]['Position']-1:tss_table.iloc[i]['Position']+20].reverse_complement()])

r_sequences35 = [
    SeqRecord(
        Seq(sequence),
        id=tss_id,
        description=""
    )
    for tss_id, sequence in sequences35
]

r_sequences10 = [
    SeqRecord(
        Seq(sequence),
        id=tss_id,
        description=""
    )
    for tss_id, sequence in sequences10
]

output_S_liv = Path(__file__).resolve().parents[2] / "data" / "S_lividans"
output_S_ros = Path(__file__).resolve().parents[2] / "data" / "S_roseum"

SeqIO.write(r_sequences35, output_S_liv / "promoters_minus35.fasta", "fasta")
SeqIO.write(r_sequences10, output_S_liv / "promoters_minus10.fasta", "fasta")

export_pr_ll(chr_plus_distances_ll, S_roseum_chromosome_seq, output_S_ros / "chr_plus_upstream.fasta")
export_pr_ll(chr_minus_distances_ll, S_roseum_chromosome_seq, output_S_ros / "chr_minus_upstream.fasta", strand='-')
export_pr_ll(pl_plus_distances_ll, S_roseum_plasmid_seq, output_S_ros / "pl_plus_upstream.fasta")
export_pr_ll(pl_minus_distances_ll, S_roseum_plasmid_seq, output_S_ros / "pl_minus_upstream.fasta", strand='-')
