"""Plot all 78 dog genes' cumulative rank distributions in one panel.

Run after rank_sum_enrichment.py: python3 src/plot_dog_gene_rank_overlay.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "notebook/results"
SERIES = (("OCD · ±1 Mb", "OCD", "1mb", "#d55e00"),
          ("OCD · ±500 kb", "OCD", "500kb", "#cc79a7"),
          ("CAD · ±1 Mb", "CAD", "1mb", "#0072b2"),
          ("CAD · ±500 kb", "CAD", "500kb", "#009e73"))


def plot():
    ranked = pd.read_csv(RESULTS / "mann_whitney_gene_ranks.csv")
    fig, ax = plt.subplots(figsize=(9, 6.5))
    max_rank = 0
    dog_ids = None
    for label, trait, window, color in SERIES:
        frame = ranked.loc[(ranked.trait == trait) & (ranked.window == window)]
        if frame.empty or frame.gene_id.duplicated().any():
            raise ValueError(f"Missing or duplicate {trait} {window} gene ranks")
        if not frame["rank"].between(1, len(frame)).all():
            raise ValueError(f"Invalid {trait} {window} gene ranks")
        dog = frame.loc[frame.dog]
        if len(dog) != 78:
            raise ValueError(f"Expected 78 dog genes in {trait} {window}; got {len(dog)}")
        ids = set(dog.gene_id)
        if dog_ids is not None and ids != dog_ids:
            raise ValueError("Dog-gene IDs differ across traits or windows")
        dog_ids = ids
        ranks = np.sort(dog["rank"].to_numpy(dtype=float))
        fraction = np.arange(1, len(ranks) + 1) / len(ranks)
        ax.step(np.r_[1, ranks, len(frame)], np.r_[0, fraction, 1],
                where="post", linewidth=2.2, color=color, label=label)
        max_rank = max(max_rank, len(frame))

    ax.plot([1, max_rank], [0, 1], color="black", linewidth=1.5,
            label="Uniform-rank reference", zorder=1)
    ax.set_xlim(1, max_rank)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([1, 5_000, 10_000, 15_000, max_rank])
    ax.set_xticklabels(["1", "5,000", "10,000", "15,000", f"{max_rank:,}"])
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_xlabel("Rank among covered protein-coding genes (strongest → weakest)")
    ax.set_ylabel("Cumulative fraction of dog genes")
    handles, labels = ax.get_legend_handles_labels()
    legend_order = [0, 2, 4, 1, 3]
    fig.legend([handles[i] for i in legend_order],
               [labels[i] for i in legend_order],
               loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.98),
               frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.83))
    fig.savefig(RESULTS / "dog_gene_rank_overlay.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    plot()
