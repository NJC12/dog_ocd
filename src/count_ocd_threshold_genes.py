"""Count confirmed dog orthologs near human OCD GWAS hits at two thresholds.

Run from the project root: python3 src/count_ocd_threshold_genes.py
"""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebook/results"
MAPPED = ROOT / "data/generated/mapped_grch38"
GENOME_WIDE = 5e-8
RELAXED = 5e-6
WINDOWS = (("1mb", 1_000_000, "min_p", "n_snps", "lead_snp", "lead_pos_grch38"),
           ("500kb", 500_000, "min_p_500kb", "n_snps_500kb", "lead_snp_500kb", "lead_pos_grch38_500kb"))


def classify(p):
    if p < GENOME_WIDE:
        return "genome_wide"
    if p < RELAXED:
        return "additional_at_5e-6"
    return "neither"


def main():
    genes = pd.read_csv(OUT / "per_gene_ocd.csv", dtype={"chrom": str})
    dog = genes.loc[genes.dog].copy()
    assert len(dog) == 78
    assert dog.dog_significant.sum() == dog.dog_suggestive.sum() == 39
    assert (dog.dog_significant ^ dog.dog_suggestive).all()
    assert dog[["n_snps", "n_snps_500kb"]].gt(0).all().all()

    rows = []
    for chrom, chrom_genes in dog.groupby("chrom"):
        variants = pd.read_parquet(MAPPED / f"chr{chrom}.parquet")
        positions = variants.pos.to_numpy()
        pvalues = variants.p.to_numpy()
        snps = variants.snp.to_numpy()
        assert np.all(positions[:-1] <= positions[1:])
        for gene in chrom_genes.itertuples(index=False):
            for window, width, p_col, n_col, snp_col, pos_col in WINDOWS:
                left = max(1, gene.tss - width)
                right = gene.tss + width
                lo = np.searchsorted(positions, left, side="left")
                hi = np.searchsorted(positions, right, side="right")
                assert hi - lo == getattr(gene, n_col)
                lead = lo + np.argmin(pvalues[lo:hi])
                assert pvalues[lead] == getattr(gene, p_col)
                assert int(positions[lead]) == int(getattr(gene, pos_col))
                assert snps[lead] == getattr(gene, snp_col)
                rows.append({
                    "window": window,
                    "dog_set": "significant" if gene.dog_significant else "suggestive",
                    "gene_id": gene.gene_id,
                    "symbol": gene.symbol,
                    "chrom": chrom,
                    "tss_grch38": gene.tss,
                    "min_p": pvalues[lead],
                    "lead_snp": snps[lead],
                    "lead_pos_grch38": int(positions[lead]),
                    "n_snps": hi - lo,
                    "threshold_category": classify(pvalues[lead]),
                })

    result = pd.DataFrame(rows).sort_values(["window", "min_p", "symbol"])
    assert len(result) == 156
    result.to_csv(OUT / "ocd_threshold_dog_genes.csv", index=False)
    summary = (
        result.groupby(["window", "dog_set", "threshold_category"])
        .size().unstack(fill_value=0)
        .reindex(columns=["genome_wide", "additional_at_5e-6", "neither"], fill_value=0)
    )
    summary["total_below_5e-6"] = summary.genome_wide + summary["additional_at_5e-6"]
    summary["genes_tested"] = summary[["genome_wide", "additional_at_5e-6", "neither"]].sum(axis=1)
    summary = summary.reset_index()
    combined = summary.groupby("window", as_index=False)[
        ["genome_wide", "additional_at_5e-6", "neither", "total_below_5e-6", "genes_tested"]
    ].sum()
    combined.insert(1, "dog_set", "combined")
    summary = pd.concat([summary, combined], ignore_index=True).sort_values(["window", "dog_set"])
    summary.to_csv(OUT / "ocd_threshold_counts.csv", index=False)
    print(summary.to_string(index=False))
    print("\nGenes below 5e-6:")
    print(result.loc[result.min_p.lt(RELAXED), ["window", "dog_set", "symbol", "min_p", "threshold_category"]].to_string(index=False))


if __name__ == "__main__":
    main()
