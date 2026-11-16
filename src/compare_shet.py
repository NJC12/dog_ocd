"""Compare human s_het estimates for all dog/human genes listed in ST6.

Run from the project root: .venv/bin/python src/compare_shet.py
"""

import gzip
import hashlib
import json
import platform
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde, linregress, mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
SHET = Path.home() / "data/constraint/shet/agarwal"
OUT = ROOT / "notebook/results/shet"
# First five aliases were verified by stable transcript ID in GENCODE v50.
# Last two: https://www.ncbi.nlm.nih.gov/gene/9289 and /gene/148137.
ALIASES = {
    "ATP5MC1": ("ATP5G1", "GENCODE transcript ENST00000393366"),
    "H4C11": ("HIST1H4J", "GENCODE transcript ENST00000355057"),
    "KATNIP": ("KIAA0556", "GENCODE transcript ENST00000261588"),
    "MAIP1": ("C2orf47", "GENCODE transcript ENST00000295079"),
    "TENT4A": ("PAPD7", "GENCODE transcript ENST00000230859"),
    "ADGRG1": ("GPR56", "NCBI Gene 9289"),
    "PROSER3": ("C19orf55", "NCBI Gene 148137"),
}


def posterior_mean(x, density):
    """E[s_het] from density in x=log10(s_het), normalized on its grid.

    Integrate 10**x * f(x) dx, not x*f(x) dx or a density in linear s.
    This is an approximation using the published, truncated KDE grid.
    """
    x, density = np.asarray(x, float), np.asarray(density, float)
    if (x.ndim != 1 or density.shape != x.shape or len(x) < 2
            or not np.isfinite(x).all() or not np.isfinite(density).all()
            or not (np.diff(x) > 0).all() or (density < 0).any()):
        raise ValueError("Invalid posterior grid or density")
    mass = trapezoid(density, x=x)
    if mass <= 0:
        raise ValueError("Posterior has no mass")
    return float(trapezoid(10.0**x * density, x=x) / mass), float(mass)


def read_posteriors(symbols, directory=SHET):
    rows, audit = [], []
    # R's syntactic column names replace the hyphen in HLA-B with a dot.
    column_to_symbol = {("HLA.B" if s == "HLA-B" else s): s for s in symbols}
    files = sorted(directory.glob("posterior_densities.*.txt.gz"))
    if len(files) != 6:
        raise ValueError("Expected six Agarwal posterior density files")
    for path in files:
        with gzip.open(path, "rt") as stream:
            header = stream.readline().rstrip().split("\t")
        selected = sorted(set(header[1:]) & set(column_to_symbol))
        if len(header) != len(set(header)):
            raise ValueError(f"Duplicate density names in {path}")
        frame = pd.read_csv(path, sep="\t", usecols=["x_axis"] + selected)
        x = frame.x_axis.to_numpy()
        if len(x) != 1024 or not np.allclose(x[[0, -1]], [-7, 0]):
            raise ValueError("Unexpected Agarwal log10(s_het) grid")
        audit.append({"file": path.name, "posterior_columns": len(header)-1,
                      "selected_columns": len(selected), "grid_points": len(x),
                      "grid_min": float(x[0]), "grid_max": float(x[-1])})
        for column in selected:
            mean, mass = posterior_mean(x, frame[column])
            rows.append({"agarwal_symbol": column_to_symbol[column],
                         "posterior_column": column, "shet_posterior_mean": mean,
                         "density_integral": mass, "posterior_file": path.name})
    result = pd.DataFrame(rows)
    if result.agarwal_symbol.duplicated().any():
        raise ValueError("Selected gene occurs in multiple posterior files")
    return result, audit


def rank_and_compare(genes):
    genes = genes.copy()
    genes["shet_rank"] = genes.shet_posterior_mean.rank(ascending=False, method="average")
    genes["log10_shet"] = np.log10(genes.shet_posterior_mean)
    dog = genes.loc[genes.group.eq("Dog"), "shet_posterior_mean"].dropna()
    human = genes.loc[genes.group.eq("Human"), "shet_posterior_mean"].dropna()
    test = mannwhitneyu(dog, human, alternative="two-sided", method="asymptotic",
                        use_continuity=True)
    # Independent rank-sum identity verifies the U statistic and its orientation.
    asc = genes.shet_posterior_mean.rank(method="average")
    u_rank = asc[genes.group.eq("Dog")].sum() - len(dog)*(len(dog)+1)/2
    if not np.isclose(test.statistic, u_rank):
        raise ValueError("Rank-sum U disagrees with scipy")
    stats = {"n_dog": len(dog), "n_human": len(human),
             "U_dog": float(test.statistic), "p_two_sided": float(test.pvalue),
             "probability_dog_greater": float(test.statistic/(len(dog)*len(human))),
             "median_shet_dog": float(dog.median()),
             "median_shet_human": float(human.median()),
             "median_rank_dog": float(genes.loc[genes.group.eq("Dog"), "shet_rank"].median()),
             "median_rank_human": float(genes.loc[genes.group.eq("Human"), "shet_rank"].median()),
             "method": "two-sided asymptotic Mann-Whitney U; tie and continuity correction"}
    return genes.sort_values(["shet_rank", "st6_symbol"], na_position="last"), stats


def plot(genes, output_dir=OUT, pdf_path=None):
    points = genes.loc[genes.shet_posterior_mean.notna() & genes.min_p_500kb.notna()]
    if not points.loc[points.group.eq("Dog"), "dog_significance"].isin(
            ["significant", "suggestive"]).all():
        raise ValueError("Missing or unexpected ST6 dog significance category")
    # Unweighted OLS in the displayed coordinates, one observation per gene.
    fit = linregress(points.neg_log10_p_500kb, points.log10_shet)
    fit_summary = {"n": len(points), "slope": float(fit.slope),
                   "intercept": float(fit.intercept), "r_squared": float(fit.rvalue**2),
                   "method": "pooled ordinary least squares with intercept in log-transformed coordinates"}
    pd.DataFrame([fit_summary]).to_csv(output_dir / "shet_pooled_regression.csv", index=False)
    print(f"Pooled OLS: slope = {fit.slope:.6f}; intercept = {fit.intercept:.6f}; n = {len(points)}")
    with plt.rc_context({"font.size": 12}):
        fig = plt.figure(figsize=(10, 8.5))
        grid = fig.add_gridspec(2, 2, width_ratios=[4, 1.15], height_ratios=[1.1, 4],
                               left=0.10, right=0.97, bottom=0.10, top=0.84,
                               hspace=0.06, wspace=0.06)
        ax = fig.add_subplot(grid[1, 0])
        ax_top = fig.add_subplot(grid[0, 0], sharex=ax)
        ax_right = fig.add_subplot(grid[1, 1], sharey=ax)
        ax_key = fig.add_subplot(grid[0, 1])
        ax_key.set_axis_off()
        series = [("Dog significant", points.group.eq("Dog") & points.dog_significance.eq("significant"), "#0072b2"),
                  ("Dog suggestive", points.group.eq("Dog") & points.dog_significance.eq("suggestive"), "#56b4e9"),
                  ("Human genes", points.group.eq("Human"), "#d55e00")]
        for label, mask, color in series:
            frame = points.loc[mask]
            ax.scatter(frame.neg_log10_p_500kb, frame.log10_shet,
                       s=45, color=color, alpha=0.8, edgecolors="white", linewidths=0.4,
                       label=f"{label} (n = {len(frame)})", zorder=3)
        for group, color in [("Dog", "#0072b2"), ("Human", "#d55e00")]:
            frame = points.loc[points.group.eq(group)]
            # Compute medians in the original units, then transform for display.
            ax.scatter(-np.log10(frame.min_p_500kb.median()),
                       np.log10(frame.shet_posterior_mean.median()),
                       marker="s", s=160, color=color, edgecolors="black",
                       linewidths=1.2, zorder=5, label=f"{group} medians")
        x_line = np.array([points.neg_log10_p_500kb.min(), points.neg_log10_p_500kb.max()])
        ax.plot(x_line, fit.intercept + fit.slope*x_line, color="#333333",
                linewidth=1.8, zorder=2, label="Pooled best fit")
        ax.text(0.03, 0.04, f"Pooled slope = {fit.slope:.3f}", transform=ax.transAxes,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9})
        ax.set_xlabel(r"$-\log_{10}(P)$, OCD ±500 kb window")
        ax.set_ylabel(r"$\log_{10}(s_{\mathrm{het}})$")
        ax.spines[["top", "right"]].set_visible(False)
        # Marginal densities are distributions of the same plotted point estimates
        # in the displayed transformed units, not mixtures of gene posteriors.
        densities = []
        for group, color in [("Dog", "#0072b2"), ("Human", "#d55e00")]:
            frame = points.loc[points.group.eq(group)]
            densities.append((group, color, gaussian_kde(frame.neg_log10_p_500kb),
                              gaussian_kde(frame.log10_shet)))
        x_pad = 3 * max(np.sqrt(kx.covariance[0, 0]) for _, _, kx, _ in densities)
        y_pad = 3 * max(np.sqrt(ky.covariance[0, 0]) for _, _, _, ky in densities)
        x_limits = (max(0, points.neg_log10_p_500kb.min() - x_pad),
                    points.neg_log10_p_500kb.max() + x_pad)
        y_limits = (points.log10_shet.min() - y_pad,
                    min(0, points.log10_shet.max() + y_pad))
        x_grid = np.linspace(*x_limits, 512)
        y_grid = np.linspace(*y_limits, 512)
        density_rows = []
        for group, color, kx, ky in densities:
            dx, dy = kx(x_grid), ky(y_grid)
            label = "Dog combined" if group == "Dog" else "Human"
            ax_top.plot(x_grid, dx, color=color, linewidth=2, label=label)
            ax_top.fill_between(x_grid, dx, color=color, alpha=0.15)
            ax_right.plot(dy, y_grid, color=color, linewidth=2)
            ax_right.fill_betweenx(y_grid, dy, color=color, alpha=0.15)
            for variable, values, density in [("neg_log10_p_500kb", x_grid, dx),
                                               ("log10_shet", y_grid, dy)]:
                density_rows.extend({"group": group, "variable": variable,
                                     "value": float(v), "density": float(d)}
                                    for v, d in zip(values, density))
        pd.DataFrame(density_rows).to_csv(output_dir / "shet_marginal_densities.csv", index=False)
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax_top.set_ylim(bottom=0)
        ax_right.set_xlim(left=0)
        ax_top.set_ylabel("Density")
        ax_right.set_xlabel("Density")
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_right.tick_params(axis="y", left=False, labelleft=False)
        for marginal in (ax_top, ax_right):
            marginal.spines[["top", "right"]].set_visible(False)
        ax_top.yaxis.set_major_locator(plt.MaxNLocator(3))
        ax_right.xaxis.set_major_locator(plt.MaxNLocator(3))
        handles, labels = ax_top.get_legend_handles_labels()
        ax_key.legend(handles, labels, title="Densities", loc="center left",
                      frameon=False, fontsize=10, title_fontsize=10)
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99),
                   ncol=3, frameon=False, fontsize=10)
        for suffix in ["png", "svg"]:
            fig.savefig(output_dir / f"shet_vs_ocd_500kb.{suffix}", dpi=200)
        if pdf_path is not None:
            fig.savefig(pdf_path, dpi=200)
        plt.close(fig)
    return len(points)


def main(pdf_path=None):
    OUT.mkdir(parents=True, exist_ok=True)
    workbook = ROOT / "data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx"
    genes = pd.read_excel(workbook, sheet_name="ST6", header=1)
    genes = genes.loc[genes.Species.isin(["Dog", "Human"]),
                      ["Gene*", "Species", "TOGA (Human:Dog)**", "dogCD significance"]].rename(
                          columns={"Gene*": "st6_symbol", "Species": "group",
                                   "TOGA (Human:Dog)**": "toga", "dogCD significance": "dog_significance"})
    genes["dog_significance"] = genes.dog_significance.str.strip().str.lower()
    if genes.st6_symbol.duplicated().any() or genes.group.value_counts().to_dict() != {"Dog": 87, "Human": 48}:
        raise ValueError("Unexpected ST6 gene membership")
    genes["human_symbol"] = genes.st6_symbol.replace({"FAM155A": "NALF1"})
    genes["agarwal_symbol"] = genes.st6_symbol.map(lambda s: ALIASES.get(s, (s, ""))[0])
    genes["shet_mapping_evidence"] = genes.st6_symbol.map(
        lambda s: ALIASES.get(s, (s, "exact ST6 symbol"))[1])
    summary_path = SHET / "elife-83172-supp2-v2.txt"
    source = pd.read_csv(summary_path, sep="\t")
    if source.Gene.duplicated().any():
        raise ValueError("Duplicate Agarwal summary symbols")
    genes = genes.merge(source[["Gene", "Ensembl_transcript_id", "group2"]].rename(
        columns={"Gene": "agarwal_symbol", "group2": "agarwal_chromosome_group"}),
        on="agarwal_symbol", how="left", validate="one_to_one")
    posteriors, density_audit = read_posteriors(genes.agarwal_symbol)
    genes = genes.merge(posteriors, on="agarwal_symbol", how="left", validate="one_to_one")
    valid = genes.shet_posterior_mean.notna()
    if not genes.loc[valid, "agarwal_chromosome_group"].eq("Autosome").all():
        raise ValueError("Expected autosomal s_het estimates")
    if not genes.loc[valid, "shet_posterior_mean"].between(1e-7, 1).all():
        raise ValueError("Posterior mean outside supplied support")
    genes["shet_status"] = np.where(valid, "available", "missing Agarwal posterior")
    windows_path = ROOT / "notebook/results/per_gene_ocd.csv"
    windows = pd.read_csv(windows_path)
    selected_windows = windows.loc[windows.symbol.isin(genes.human_symbol),
                                  ["symbol", "gene_id", "min_p_500kb", "n_snps_500kb", "dog"]]
    genes = genes.merge(selected_windows.rename(columns={"symbol": "human_symbol", "dog": "in_previous_78_dog_set"}),
                        on="human_symbol", how="left", validate="one_to_one")
    if genes.gene_id.isna().any() or not genes.min_p_500kb.between(0, 1, inclusive="right").all():
        raise ValueError("Missing human mapping or invalid/nonpositive OCD window p")
    genes["neg_log10_p_500kb"] = -np.log10(genes.min_p_500kb)
    genes, stats = rank_and_compare(genes)
    genes.to_csv(OUT / "shet_gene_ranks.csv", index=False)
    pd.DataFrame([stats]).to_csv(OUT / "shet_mann_whitney.csv", index=False)
    stats["n_plotted"] = plot(genes, pdf_path=pdf_path)
    stats["missing_shet_genes"] = genes.loc[genes.shet_posterior_mean.isna(), "st6_symbol"].tolist()
    input_paths = [workbook, windows_path, summary_path] + sorted(SHET.glob("posterior_densities.*.txt.gz"))
    audit = {"results": stats, "density_files": density_audit,
             "density_integral_range": [float(genes.density_integral.min()), float(genes.density_integral.max())],
             "environment": {"python": platform.python_version(), "platform": platform.platform(),
                             "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
             "sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths}}
    (OUT / "shet_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print("All ST6 genes: 87 Dog + 48 Human; human s_het estimates for both groups.")
    print("Rank 1 = highest posterior mean s_het; average ranks for ties.")
    print(pd.DataFrame([stats]).to_string(index=False))
    print("\nRanked genes:")
    print(genes[["shet_rank", "st6_symbol", "group", "shet_posterior_mean", "log10_shet", "min_p_500kb"]].to_string(index=False))
    print(f"\nPlot: {OUT / 'shet_vs_ocd_500kb.png'}")
    return genes, stats


if __name__ == "__main__":
    main()
