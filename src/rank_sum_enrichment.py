"""Rank per-gene OCD and CAD signals and test dog-gene enrichment.

Run after compare_ocd.py and compare_cad.py: python3 src/rank_sum_enrichment.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'notebook/results'
WINDOWS = (('1mb', '', 1_000_000), ('500kb', '_500kb', 500_000))
DOG_SETS = (('significant', 'dog_significant'),
            ('suggestive', 'dog_suggestive'),
            ('combined', 'dog'))


def rank_and_test(frame, trait, window, suffix, width):
    """Rank covered genes and test three dog sets against the same non-dog genes."""
    required = {'gene_id', 'symbol', 'chrom', 'dog', 'dog_significant',
                'dog_suggestive', f'min_p{suffix}', f'n_snps{suffix}'}
    if missing := required.difference(frame):
        raise ValueError(f'{trait} is missing columns: {sorted(missing)}')
    if frame.gene_id.duplicated().any():
        raise ValueError(f'{trait} has duplicate gene IDs')
    if not frame.dog.eq(frame.dog_significant | frame.dog_suggestive).all():
        raise ValueError(f'{trait} dog labels do not match significant/suggestive union')
    if (frame.dog_significant & frame.dog_suggestive).any():
        raise ValueError(f'{trait} significant and suggestive dog sets overlap')

    p_col, count_col = f'min_p{suffix}', f'n_snps{suffix}'
    covered = frame.loc[frame[count_col].gt(0),
                        ['gene_id', 'symbol', 'chrom', 'dog', 'dog_significant',
                         'dog_suggestive', p_col, count_col]].copy()
    p = covered[p_col].to_numpy(dtype=float)
    if not np.isfinite(p).all() or ((p <= 0) | (p > 1)).any():
        raise ValueError(f'{trait} {window} has invalid covered-gene p-values')
    covered = covered.rename(columns={p_col: 'min_p', count_col: 'n_snps'})
    covered.insert(0, 'trait', trait)
    covered.insert(1, 'window', window)
    covered.insert(2, 'width_each_side_bp', width)
    covered['max_neg_log10_p'] = -np.log10(p)
    covered['rank'] = covered.max_neg_log10_p.rank(
        ascending=False, method='average')
    covered = covered.sort_values(['rank', 'gene_id'], kind='mergesort')

    other = covered.loc[~covered.dog, 'max_neg_log10_p'].to_numpy()
    summaries = []
    for set_name, mask_col in DOG_SETS:
        dog = covered.loc[covered[mask_col], 'max_neg_log10_p'].to_numpy()
        if len(dog) == 0 or len(other) == 0:
            raise ValueError(f'{trait} {window} {set_name} needs both dog and other covered genes')
        test = mannwhitneyu(dog, other, alternative='greater', method='asymptotic')
        summaries.append(dict(
            trait=trait, window=window, dog_set=set_name,
            width_each_side_bp=width, n_dog=len(dog), n_other=len(other),
            median_dog_min_p=float(np.median(covered.loc[covered[mask_col], 'min_p'])),
            median_other_min_p=float(np.median(covered.loc[~covered.dog, 'min_p'])),
            median_dog_rank=float(covered.loc[covered[mask_col], 'rank'].median()),
            median_other_rank=float(covered.loc[~covered.dog, 'rank'].median()),
            u_statistic=float(test.statistic),
            probability_dog_stronger=float(test.statistic / (len(dog) * len(other))),
            one_sided_p=float(test.pvalue)))
    return covered, summaries


def main():
    ranks, summaries = [], []
    for trait in ('OCD', 'CAD'):
        frame = pd.read_csv(OUT / f'per_gene_{trait.lower()}.csv', dtype={'chrom': str})
        for window, suffix, width in WINDOWS:
            ranked, window_summaries = rank_and_test(frame, trait, window, suffix, width)
            ranks.append(ranked)
            summaries.extend(window_summaries)
    results = pd.DataFrame(summaries)
    results['holm_adjusted_p'] = multipletests(results.one_sided_p, method='holm')[1]
    pd.concat(ranks, ignore_index=True).to_csv(OUT / 'mann_whitney_gene_ranks.csv', index=False)
    results.to_csv(OUT / 'mann_whitney_results.csv', index=False)
    print(results.to_string(index=False))


if __name__ == '__main__':
    main()
