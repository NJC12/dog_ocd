"""Repeat dog-gene enrichment using gene-level OCD pLoF SKAT-O results."""

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/downloaded/burden/phenotype-gene-burden_obsessive_compulsive_disorder_meta_pLoF_maf-0.001_2026_10_15_15_55_32.csv'
OUT = ROOT / 'notebook/results/burden'
GROUPS = (('combined', 'dog', 'All dog genes', '#4d4d4d'),
          ('significant', 'dog_significant', 'Dog GWAS significant', '#d55e00'),
          ('suggestive', 'dog_suggestive', 'Dog GWAS suggestive', '#0072b2'))


def analyze(genes, burden):
    if genes.gene_id.duplicated().any() or burden['Gene ID'].duplicated().any():
        raise ValueError('Duplicate gene IDs')
    p = burden['P SKAT-O'].to_numpy(dtype=float)
    if not np.isfinite(p).all() or not ((p > 0) & (p <= 1)).all():
        raise ValueError('Invalid SKAT-O p-values')
    cols = ['gene_id', 'symbol', 'chrom', 'dog', 'dog_significant', 'dog_suggestive']
    ranked = genes[cols].merge(burden, left_on='gene_id', right_on='Gene ID', validate='one_to_one')
    assert ranked.dog.eq(ranked.dog_significant | ranked.dog_suggestive).all()
    assert not (ranked.dog_significant & ranked.dog_suggestive).any()
    ranked['rank'] = ranked['P SKAT-O'].rank(method='average')
    ranked = ranked.sort_values(['rank', 'gene_id'])
    other = ranked.loc[~ranked.dog, 'P SKAT-O']
    rows = []
    for name, col, _, _ in GROUPS:
        dog = ranked.loc[ranked[col], 'P SKAT-O']
        # Higher -P is stronger: matches the established enrichment U orientation.
        test = mannwhitneyu(-dog, -other, alternative='greater', method='asymptotic')
        two = mannwhitneyu(-dog, -other, alternative='two-sided', method='asymptotic')
        pooled = pd.Series(np.r_[-dog.to_numpy(), -other.to_numpy()]).rank()
        independent_u = pooled.iloc[:len(dog)].sum() - len(dog) * (len(dog) + 1) / 2
        assert np.isclose(independent_u, test.statistic)
        rows.append(dict(dog_set=name, n_dog=len(dog), n_other=len(other),
                         median_dog_rank=ranked.loc[ranked[col], 'rank'].median(),
                         median_dog_p=dog.median(), median_other_p=other.median(),
                         u_statistic=test.statistic,
                         probability_dog_stronger=test.statistic / (len(dog) * len(other)),
                         one_sided_p=test.pvalue, two_sided_p=two.pvalue))
    results = pd.DataFrame(rows)
    results['holm_adjusted_p'] = multipletests(results.one_sided_p, method='holm')[1]
    return ranked, results


def plot(ranked):
    for kind in ('rank', 'p'):
        fig, ax = plt.subplots(figsize=(9, 6.5))
        lo, hi = (1, len(ranked)) if kind == 'rank' else (np.floor(np.log10(ranked['P SKAT-O'].min())), 0)
        for _, col, label, color in GROUPS:
            group = ranked.loc[ranked[col]]
            values = np.sort(group['rank'] if kind == 'rank' else np.log10(group['P SKAT-O']))
            y = np.arange(1, len(values) + 1) / len(values)
            ax.step(np.r_[lo, values, hi], np.r_[0, y, 1], where='post',
                    color=color, linewidth=2.2, label=f'{label} (n={len(values)})')
        if kind == 'rank':
            ax.plot([lo, hi], [0, 1], color='black', linewidth=1.5, label='Uniform-rank reference', zorder=1)
        ax.set(xlim=(lo, hi), ylim=(0, 1.02), ylabel='Cumulative fraction of dog genes',
               xlabel='Rank among tested autosomal protein-coding genes (strongest → weakest)'
               if kind == 'rank' else 'log₁₀(SKAT-O P cutoff)')
        ax.set_yticks(np.linspace(0, 1, 5))
        fig.legend(*ax.get_legend_handles_labels(), loc='upper center', ncol=2,
                   bbox_to_anchor=(0.5, 0.98), frameon=False)
        fig.tight_layout(rect=(0, 0, 1, 0.85))
        for ext in ('png', 'svg'):
            fig.savefig(OUT / f'dog_gene_{kind}_cdfs.{ext}', dpi=200)
        plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / 'notebook/results/per_gene_ocd.csv', dtype={'chrom': str})
    burden = pd.read_csv(SOURCE)
    ranked, results = analyze(genes, burden)
    missing = genes.loc[genes.dog & ~genes.gene_id.isin(burden['Gene ID'])]
    missing.to_csv(OUT / 'missing_dog_genes.csv', index=False)
    ranked.to_csv(OUT / 'gene_ranks.csv', index=False)
    results.to_csv(OUT / 'mann_whitney_results.csv', index=False)
    audit = dict(source=SOURCE.name, sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                 input_rows=len(burden), analyzed_genes=len(ranked),
                 excluded_outside_original_background=len(burden)-len(ranked),
                 missing_gene_symbols=int(burden.Gene.isna().sum()),
                 missing_dog_genes=missing.symbol.tolist(),
                 min_p=float(burden['P SKAT-O'].min()), max_p=float(burden['P SKAT-O'].max()))
    (OUT / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    plot(ranked)
    report = results.to_string(index=False)
    (OUT / 'analysis_output.txt').write_text(report+'\n')
    print(report)
    return ranked, results


if __name__ == '__main__':
    main()
