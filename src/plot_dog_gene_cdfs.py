"""Plot cumulative human GWAS p-value distributions for confirmed dog genes.

Run from the project root: python3 src/plot_dog_gene_cdfs.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "notebook/results"
WINDOWS = (("±1 Mb", "min_p", "n_snps"),
           ("±500 kb", "min_p_500kb", "n_snps_500kb"))
GROUPS = (("All dog genes", "dog", "#4d4d4d"),
          ("Dog GWAS significant", "dog_significant", "#d55e00"),
          ("Dog GWAS suggestive", "dog_suggestive", "#0072b2"))


def plot():
    data = {trait: pd.read_csv(RESULTS / f"per_gene_{trait.lower()}.csv")
            for trait in ("OCD", "CAD")}
    for frame in data.values():
        assert frame.dog.sum() == 78
        assert frame.dog_significant.sum() == frame.dog_suggestive.sum() == 39
        assert (frame.loc[frame.dog, "dog_significant"] ^
                frame.loc[frame.dog, "dog_suggestive"]).all()

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    x_min = -14
    for row, (trait, frame) in enumerate(data.items()):
        for col, (window, p_col, count_col) in enumerate(WINDOWS):
            ax = axes[row, col]
            for label, mask_col, color in GROUPS:
                group = frame.loc[frame[mask_col] & frame[count_col].gt(0), p_col]
                values = group.to_numpy(dtype=float)
                if not len(values) or not np.isfinite(values).all() or not ((values > 0) & (values <= 1)).all():
                    raise ValueError(f"Invalid or missing {trait} {window} p-values for {label}")
                # At cutoff 10**x, the curve equals the fraction with a
                # window minimum p-value no greater than that cutoff.
                logp = np.sort(np.log10(values))
                y = np.arange(1, len(logp) + 1) / len(logp)
                ax.step(np.r_[x_min, logp, 0], np.r_[0, y, 1],
                        where="post", color=color, linewidth=2,
                        label=f"{label} (n={len(logp)})")
            ax.set_title(f"{trait} · {window}")
            ax.set_xlim(x_min, 0)
            ax.set_ylim(0, 1.02)
            ax.set_xticks(np.arange(-14, 1, 2))
            ax.set_yticks(np.linspace(0, 1, 5))
            ax.tick_params(axis="x", labelbottom=True)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.98), frameon=False)
    fig.supxlabel("log₁₀(P cutoff)")
    fig.supylabel("Cumulative fraction of dog genes")
    fig.tight_layout(rect=(0.03, 0.03, 1, 0.92))
    fig.savefig(RESULTS / "dog_gene_p_cdfs.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    plot()
