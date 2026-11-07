import json


def operon_search(promoters: list, pr_groups: list) -> list:
    a, b = 0, 0
    while not (a == len(promoters) - 1 and b == len(pr_groups)-1):
        if (pr_groups[b][0][1] - promoters[a][1] >= 6 and pr_groups[b][0][1] - promoters[a][1] <= 20):
            yield (pr_groups[b], promoters[a])
        if (pr_groups[b][0][1] - promoters[a][1] < 6 or a == len(promoters) - 1) and b != len(pr_groups) - 1:
            b+=1
        else:
            a+=1

with open("jsons/chr_plus_distances.json", "r") as f:
    chr_plus_d = json.load(f)

with open("jsons/pred_pr_chr_plus_chain.json", "r") as f:
    chr_plus_pr = json.load(f)

with open("test.json", "w") as f:
    json.dump(list(operon_search(chr_plus_pr, chr_plus_d)), f)

print(len(list(operon_search(chr_plus_pr, chr_plus_d))))