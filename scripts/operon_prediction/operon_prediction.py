import json


def operon_search(promoters: list, pr_groups: list) -> list:
    a, b = 0, 0
    while not (a == len(promoters) - 1 and b == len(pr_groups)-1):
        if (pr_groups[b][0][1] - promoters[a][1] >= 25 and pr_groups[b][0][1] - promoters[a][1] <= 41):
            yield (pr_groups[b], promoters[a])
        if (pr_groups[b][0][1] - promoters[a][1] < 25 or a == len(promoters) - 1) and b != len(pr_groups) - 1:
            b+=1
        else:
            a+=1
    if (pr_groups[b][0][1] - promoters[a][1] >= 25 and pr_groups[b][0][1] - promoters[a][1] <= 41):
        yield (pr_groups[b], promoters[a])
    

def operons_on_seq(dist_file_name: str, pr_file_name: str, output_name: str):
    with open(dist_file_name, "r") as f:
        chr_plus_d = json.load(f)

    with open(pr_file_name, "r") as f:
        chr_plus_pr = json.load(f)

    with open(output_name, "w") as f:
        json.dump(list(operon_search(chr_plus_pr, chr_plus_d)), f)
        print(len(list(operon_search(chr_plus_pr, chr_plus_d))))


operons_on_seq("jsons/chr_plus_distances.json", "jsons/pred_pr_chr_plus_chain.json", "jsons/res/chr_plus.json")
operons_on_seq("jsons/chr_minus_distances.json", "jsons/pred_pr_chr_minus_chain.json", "jsons/res/chr_minus.json")
operons_on_seq("jsons/pl_plus_distances.json", "jsons/pred_pr_pl_plus_chain.json", "jsons/res/pl_plus.json")
operons_on_seq("jsons/pl_minus_distances.json", "jsons/pred_pr_pl_minus_chain.json", "jsons/res/pl_minus.json")

# print(len(list(operon_search(chr_plus_pr, chr_plus_d))))
# results = list(operon_search(chr_plus_pr, chr_plus_d))

# unique_groups = {
#     (group[0][1], group[-1][2])
#     for group, promoter in results
# }

# print("Всего совпадений:", len(results))
# print("Уникальных групп:", len(unique_groups))

# count = []
# u = 0
# for idx, i in enumerate(results):
#     a = i[0][0][1]
#     b = i[0][-1][2]
#     while u < len(chr_plus_pr) and chr_plus_pr[u][1] < b:
#         if (chr_plus_pr[u][1] <= b or chr_plus_pr[u][0] >= a) and idx not in count:
#             count.append(idx)
#         u+=1
#     # if (chr_plus_pr[u][1] <= b or chr_plus_pr[u][0] >= a) and idx not in count:
#     #     count.append(idx)
#     if u >= len(chr_plus_pr):
#         break
# print(len(count))