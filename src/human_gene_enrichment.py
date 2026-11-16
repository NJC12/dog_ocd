"""Enrichment of the ST6 Human set using the established GWAS gene windows."""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def read_human_gene_set(workbook, genes):
    """Map all 48 ST6 Human symbols exactly to the annotation used for scoring."""
    st6 = pd.read_excel(workbook, sheet_name='ST6', header=1)
    human = st6.loc[st6.Species.eq('Human'), ['Gene*']].rename(
        columns={'Gene*': 'st6_symbol'})
    if len(human) != 48 or human.st6_symbol.isna().any() or human.st6_symbol.duplicated().any():
        raise ValueError('Expected 48 unique Human genes in ST6')
    mapping = genes.loc[genes.symbol.isin(human.st6_symbol), ['gene_id', 'symbol']]
    human = human.merge(mapping, left_on='st6_symbol', right_on='symbol',
                        how='left', validate='one_to_one')
    if human.gene_id.isna().any() or human.gene_id.duplicated().any():
        raise ValueError('Missing or ambiguous ST6 Human gene mapping')
    return human


def human_mann_whitney(ranked, human_ids):
    """Test one human set against its complement; Holm spans all supplied groups."""
    human_ids = set(human_ids)
    if not human_ids:
        raise ValueError('The Human gene set is empty')
    rows = []
    for (trait, window), frame in ranked.groupby(['trait', 'window'], sort=False):
        if frame.gene_id.duplicated().any():
            raise ValueError(f'Duplicate gene IDs in {trait} {window}')
        selected = frame.gene_id.isin(human_ids)
        if set(frame.loc[selected, 'gene_id']) != human_ids:
            raise ValueError(f'Not all Human genes are covered in {trait} {window}')
        p = frame.min_p.to_numpy(dtype=float)
        if not np.isfinite(p).all() or not ((p > 0) & (p <= 1)).all():
            raise ValueError(f'Invalid window p-values in {trait} {window}')
        human, other = frame.loc[selected], frame.loc[~selected]
        if other.empty:
            raise ValueError('The non-Human-set background is empty')
        test = mannwhitneyu(-np.log10(human.min_p), -np.log10(other.min_p),
                           alternative='greater', method='asymptotic')
        rows.append({
            'trait': trait, 'window': window, 'gene_set': 'ST6 Human',
            'width_each_side_bp': int(frame.width_each_side_bp.iloc[0]),
            'n_human': len(human), 'n_other': len(other),
            'median_human_min_p': float(human.min_p.median()),
            'median_other_min_p': float(other.min_p.median()),
            'median_human_rank': float(human['rank'].median()),
            'median_other_rank': float(other['rank'].median()),
            'u_statistic': float(test.statistic),
            'probability_human_stronger': float(test.statistic / (len(human) * len(other))),
            'one_sided_p': float(test.pvalue),
        })
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError('No trait/window groups supplied')
    result['holm_adjusted_p'] = multipletests(result.one_sided_p, method='holm')[1]
    result['holm_family_tests'] = len(result)
    return result


def circular_shift_gene_set(frame, selected_ids, n=10_000, p_col='min_p',
                            count_col='n_snps', seed=20261013):
    """Repeat the original circular shifts for one gene set.

    Draw offsets in the original replicate/chromosome order. Indexed lookup of
    shifted selected positions avoids repeatedly rolling full boolean arrays;
    positive offsets have exactly the same orientation as numpy.roll(labels).
    """
    if n < 1 or frame.gene_id.duplicated().any():
        raise ValueError('Need positive shift count and unique gene IDs')
    covered = frame.loc[frame[count_col].gt(0)].copy()
    selected_ids = set(selected_ids)
    covered['human'] = covered.gene_id.isin(selected_ids)
    if not selected_ids or set(covered.loc[covered.human, 'gene_id']) != selected_ids:
        raise ValueError('All Human genes must have covered windows')
    p = covered[p_col].to_numpy(dtype=float)
    if not np.isfinite(p).all() or not ((p > 0) & (p <= 1)).all():
        raise ValueError('Invalid covered-gene p-values')
    covered['snp_stratum'] = covered.groupby('chrom')[count_col].transform(
        lambda x: pd.qcut(x.rank(method='first'), min(5, len(x)), labels=False))
    covered['signal_rank'] = covered.groupby(['chrom', 'snp_stratum'])[p_col].transform(
        lambda x: x.rank(ascending=True, method='average', pct=True))
    chroms = [group.sort_values('tss') for _, group in covered.groupby('chrom')]
    rng = np.random.default_rng(seed)
    offsets = np.array([[rng.integers(len(group)) for group in chroms]
                        for _ in range(n)])
    shifted_sums = np.zeros(n)
    for j, group in enumerate(chroms):
        positions = np.flatnonzero(group.human.to_numpy())
        if len(positions):
            scores = group.signal_rank.to_numpy()
            shifted = (positions[None, :] + offsets[:, j, None]) % len(group)
            shifted_sums += scores[shifted].sum(axis=1)
    null = shifted_sums / len(selected_ids)
    observed = float(covered.loc[covered.human, 'signal_rank'].mean())
    # Treat roundoff in an equal sum as equality, as required by the <= test.
    null_le_observed = (null < observed) | np.isclose(null, observed, rtol=0, atol=1e-14)
    stats = {
        'n': len(selected_ids), 'observed_mean_rank': observed,
        'null_mean': float(null.mean()),
        'empirical_one_sided_p': float((1 + null_le_observed.sum()) / (n + 1)),
        'permutations': n, 'seed': seed,
    }
    return covered, stats
