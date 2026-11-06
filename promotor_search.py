import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def find_seq_fast(seq: Seq, target: str, strain: str, max_mismatches=2):
    arr = np.frombuffer(
        str(seq).encode("ascii"),
        dtype=np.uint8
    )

    motif = np.frombuffer(
        target.encode("ascii"),
        dtype=np.uint8
    )

    k = len(motif)
    errors = np.zeros(len(arr) - k + 1, dtype=np.uint8)

    for i in range(k):
        errors += (arr[i:i + len(errors)] != motif[i])

    positions = np.flatnonzero(errors <= max_mismatches)

    return [
        (int(pos), int(pos + k), strain, int(errors[pos]))
        for pos in positions
    ]


def both_seq(list_1: list, list_2: list):
    a, b = 0, 0
    while not (a == len(list_1) - 1 and b == len(list_2) - 1):
        if (list_2[b][0] - list_1[a][0] <= 30) and (list_2[b][0] - list_1[a][0] >= 10):
            yield (list_1[a][0], list_2[b][0])
        if (list_2[b][0] - list_1[a][0] > 30 and a != len(list_1)-1) or b == len(list_2)-1:
            a+=1
        else:
            b+=1


sequences = []

for record in SeqIO.parse("genomic.fna", "fasta"):
    print(record.id)
    print(len(record.seq))
    print(record.description)
    sequences.append(record.seq)

chromosome_seq = sequences[0].upper()
plasmid_seq = sequences[1].upper()

start = 0
target = 'TATAAT'
Pribnow_box_res = list(find_seq_fast(seq=chromosome_seq,target='TATAAT',strain='+'))
el_35 = list(find_seq_fast(seq=chromosome_seq, target='TTGACA',strain='+'))


print(len(list(both_seq(el_35, Pribnow_box_res))))