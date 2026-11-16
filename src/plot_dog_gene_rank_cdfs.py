"""Plot dog-gene cumulative distributions over genome-wide gene-signal ranks.

Run after rank_sum_enrichment.py: python3 src/plot_dog_gene_rank_cdfs.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "notebook/results"
WINDOWS = (("1mb", "±1 Mb"), ("500kb", "±500 kb"))
GROUPS = (("All dog genes", "dog", "#4d4d4d", 78),
          ("Dog GWAS significant", "dog_significant", "#d55e00", 39),
          ("Dog GWAS suggestive", "dog_suggestive", "#0072b2", 39))


def plot():
    ranked = pd.read_csv(RESULTS / "mann_whitney_gene_ranks.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    max_rank = 0
    for row, trait in enumerate(("OCD", "CAD")):
        for col, (window, title) in enumerate(WINDOWS):
            frame = ranked.loc[(ranked.trait == trait) & (ranked.window == window)]
            if frame.empty or frame.gene_id.duplicated().any():
                raise ValueError(f"Missing or duplicate {trait} {window} gene ranks")
            if not frame["rank"].between(1, len(frame)).all():
                raise ValueError(f"Invalid {trait} {window} gene ranks")
            if not frame.dog.eq(frame.dog_significant | frame.dog_suggestive).all():
                raise ValueError(f"Inconsistent {trait} {window} dog labels")
            if (frame.dog_significant & frame.dog_suggestive).any():
                raise ValueError(f"Overlapping {trait} {window} dog sets")
            max_rank = max(max_rank, len(frame))
            ax = axes[row, col]
            for label, mask_col, color, expected_n in GROUPS:
                ranks = np.sort(frame.loc[frame[mask_col], "rank"].to_numpy(dtype=float))
                if len(ranks) != expected_n:
                    raise ValueError(f"Expected {expected_n} {label} in {trait} {window}; got {len(ranks)}")
                fraction = np.arange(1, len(ranks) + 1) / len(ranks)
                ax.step(np.r_[1, ranks, len(frame)], np.r_[0, fraction, 1],
                        where="post", color=color, linewidth=2,
                        label=f"{label} (n={len(ranks)})")
            ax.set_title(f"{trait} · {title}")
            ax.set_ylim(0, 1.02)
            ax.set_yticks(np.linspace(0, 1, 5))
            ax.tick_params(axis="x", labelbottom=True)

    for ax in axes.flat:
        ax.set_xlim(1, max_rank)
        ax.set_xticks([1, 5_000, 10_000, 15_000, max_rank])
        ax.set_xticklabels(["1", "5,000", "10,000", "15,000", f"{max_rank:,}"])
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.98), frameon=False)
    fig.supxlabel("Rank among covered protein-coding genes (strongest → weakest)")
    fig.supylabel("Cumulative fraction of dog genes")
    fig.tight_layout(rect=(0.03, 0.03, 1, 0.92))
    fig.savefig(RESULTS / "dog_gene_rank_cdfs.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    plot()
