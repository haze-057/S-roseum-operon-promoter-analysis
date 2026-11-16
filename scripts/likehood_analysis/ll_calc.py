import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO


def intercds_distances(cdss: pd.DataFrame, pairs: pd.DataFrame,) -> tuple[np.ndarray, np.ndarray]:
    within_operon = np.array([])
    not_within_operon = np.array([])
    for i in range(pairs.shape[0]):
        if pairs.iloc[i]['pair_type'] == 'KOP':
            within_operon = np.append(within_operon, cdss[cdss['gene'] == pairs.iloc[i]['gene_2']].iloc[0]['start'] - cdss[cdss['gene'] == pairs.iloc[i]['gene_1']].iloc[0]['end'] )
        else:
            not_within_operon = np.append(not_within_operon, cdss[cdss['gene'] == pairs.iloc[i]['gene_2']].iloc[0]['start'] - cdss[cdss['gene'] == pairs.iloc[i]['gene_1']].iloc[0]['end'])
    return within_operon, not_within_operon        

def likehood_calculation(within_op: np.ndarray, not_within_op: np.ndarray, dist: int, n_bins: int, a: float = 0.5, w_bins: int = 10) -> float:
    # return np.log((np.mean(within_op == dist) + 0.5) / (np.mean(not_within_op == dist) + 0.5 * 10))
    return (np.log(((np.count_nonzero((within_op >= dist) & (within_op < dist + w_bins)) + a) / (len(within_op) + a * n_bins)) / ((np.count_nonzero((not_within_op >= dist) & (not_within_op < dist + w_bins)) + a) / (len(not_within_op) + a * n_bins))))

def plot_calc(within_op: np.ndarray, not_within_op: np.ndarray, w_bins: int = 10) -> np.ndarray:
    ans = np.array([])
    min_dist = min(int(within_op.min()), int(not_within_op.min()))
    max_dist = max(int(within_op.max()), int(not_within_op.max()))
    for i in range(min_dist, max_dist - 1 + w_bins, w_bins):
        ans = np.append(ans, likehood_calculation(within_op, not_within_op, i, n_bins=((max_dist - min_dist - 1) // w_bins + 1), w_bins=w_bins))
    return ans, min_dist, max_dist

def get_ll(dist: int, res: np.ndarray, min_dist: int, w_bins: int = 10):
    index = (dist - min_dist) // w_bins
    if index < 0 or index >= len(res):
        return np.nan
    return res[int(index)]

def figure_1(within_operon: np.ndarray, not_within_operon: np.ndarray):
    plt.figure(figsize=(12, 5))
    bins = np.histogram_bin_edges(np.concatenate((within_operon, not_within_operon)), bins=50)
    centers = (bins[:-1] + bins[1:]) / 2
    kernel = np.exp(-0.5 * (np.arange(-6, 7) / 2) ** 2)
    kernel /= kernel.sum()

    for distances, label in ((within_operon, "Within operon"), (not_within_operon, "Outside operon")):
        density, _ = np.histogram(distances, bins=bins, density=True)
        smooth_density = np.convolve(np.pad(density, (6, 6), mode="reflect"), kernel, mode="valid")
        plt.plot(centers, smooth_density, label=label)
        plt.fill_between(centers, smooth_density, alpha=0.2)

    plt.xlabel("Inter-CDS distance (bp)")
    plt.ylabel("Probability density")
    plt.legend()

    plt.tight_layout()
    plt.savefig(results_dir / "Figure_1.png", dpi=300)
    plt.close()

def figure_2(within_operon: np.ndarray, not_within_operon: np.ndarray):
    plt.figure(figsize=(12, 5))
    plot_res = res[:252]
    distances = min_dist + np.arange(plot_res.size) * 10

    plt.bar(
        distances,
        plot_res,
        width=10,
        align="edge",
        edgecolor="black",
        linewidth=0.3
    )

    plt.axhline(0, linestyle="--", color="gray")

    plt.xlabel("Inter-CDS distance (bp)")
    plt.ylabel("Log-likelihood ratio")
    plt.tight_layout()
    plt.savefig(results_dir / "Figure_2.png", dpi=300)
    plt.close()

def figure_3(within_operon: np.ndarray, not_within_operon: np.ndarray):
    thresholds = np.linspace(-3, 3, 100)

    kop_scores = np.array([get_ll(dist, res, min_dist, w_bins=10) for dist in within_operon])
    nop_scores = np.array([get_ll(dist, res, min_dist, w_bins=10) for dist in not_within_operon])

    kop_ratio = np.array([
        np.mean(kop_scores >= t)
        for t in thresholds
    ])

    nop_ratio = np.array([
        np.mean(nop_scores >= t)
        for t in thresholds
    ])

    j_scores = kop_ratio - nop_ratio

    best_index = np.argmax(j_scores)
    best_threshold = thresholds[best_index]

    print("Best threshold:", best_threshold)
    print("Youden's J:", j_scores[best_index])
    
    plt.figure(figsize=(10, 6))

    plt.plot(thresholds, kop_ratio, label="KOP retained")
    plt.plot(thresholds, nop_ratio, label="NOP incorrectly retained")
    plt.plot(thresholds, j_scores, label="Youden's J", linestyle="--")

    plt.axvline(
        best_threshold,
        linestyle=":",
        color="gray",
        label=f"Optimal threshold: {best_threshold:.2f}"
    )

    plt.xlabel("LL(dist) threshold")
    plt.ylabel("Fraction / Youden's J")
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(results_dir / "Figure_3.png", dpi=300)
    plt.close()

project_root = Path(__file__).resolve().parents[2]
data_dir = project_root / "data" / "S_coelicolor"
results_dir = project_root / "results" / "operon_prediction_ll"
results_dir.mkdir(parents=True, exist_ok=True)

record = SeqIO.read(data_dir / "sequence.gb", "genbank")

genes = []

for feature in record.features:
    if feature.type == "CDS":
        genes.append({
            "gene": feature.qualifiers.get("gene", [None])[0],
            "start": int(feature.location.start),
            "end": int(feature.location.end),
            "strand": "+" if feature.location.strand == 1 else "-",
            "product": feature.qualifiers.get("product", [None])[0]
        })

gb_table = pd.DataFrame(genes).sort_values("start")

pr_pairs = pd.read_csv(data_dir / "gene_pairs.CSV")

within_operon, not_within_operon = intercds_distances(gb_table, pr_pairs)

res, min_dist, max_dist = plot_calc(within_op=within_operon,not_within_op=not_within_operon)

with open(results_dir / "likehood.json", "w") as f:
    json.dump((min_dist, max_dist, res.tolist()), f)

figure_1(within_operon, not_within_operon)

figure_2(within_operon, not_within_operon)

figure_3(within_operon, not_within_operon)

