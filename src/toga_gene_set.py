"""Select dog genes using only the supplementary workbook's TOGA labels."""

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def adjust_ocd_tests(results):
    """Holm correction for the six OCD dog-set tests; controls stay unadjusted."""
    results = results.copy()
    primary = results.trait.eq('OCD')
    expected = {(window, group) for window in ['1mb', '500kb']
                for group in ['significant', 'suggestive', 'combined']}
    family = results.loc[primary]
    if len(family) != 6 or set(zip(family.window, family.dog_set)) != expected:
        raise ValueError('Expected exactly six unique OCD window-by-dog-set tests')
    if not family.one_sided_p.between(0, 1).all():
        raise ValueError('Invalid OCD p-values')
    results['holm_adjusted_p'] = np.nan
    results.loc[primary, 'holm_adjusted_p'] = multipletests(family.one_sided_p, method='holm')[1]
    results['holm_family_tests'] = np.where(primary, 6, 0)
    return results


def read_toga_gene_set(workbook, genes):
    """Map ST6 one-to-one Dog symbols to the human scoring annotation.

    FAM155A is the former name of NALF1. No Ensembl homology decisions
    participate in selection; stable IDs identify the human scoring windows.
    """
    st6 = pd.read_excel(workbook, sheet_name='ST6', header=1)
    dog = st6.loc[st6.Species.eq('Dog')].copy()
    dog['TOGA'] = dog['TOGA (Human:Dog)**'].astype(str).str.strip()
    selected = dog.loc[dog.TOGA.eq('1:1')].copy()
    selected['human_symbol'] = selected['Gene*'].replace({'FAM155A': 'NALF1'})
    mapping = genes.loc[genes.symbol.isin(selected.human_symbol),
                        ['gene_id', 'symbol']].rename(
                            columns={'gene_id': 'human_gene_id', 'symbol': 'human_symbol'})
    selected = selected.merge(mapping, on='human_symbol', how='left', validate='one_to_one')
    if selected.human_gene_id.isna().any() or selected.human_gene_id.duplicated().any():
        raise ValueError('Missing or ambiguous TOGA human gene mapping')
    if not selected['dogCD significance'].isin(['significant', 'suggestive']).all():
        raise ValueError('Unexpected dog significance label')
    return dog, selected


def apply_dog_labels(frame, selected):
    """Replace membership labels while preserving all GWAS window scores."""
    frame = frame.copy()
    if not set(selected.human_gene_id).issubset(set(frame.gene_id)):
        raise ValueError('TOGA genes missing from window scores')
    frame['dog'] = frame.gene_id.isin(selected.human_gene_id)
    for name in ['significant', 'suggestive']:
        ids = selected.loc[selected['dogCD significance'].eq(name), 'human_gene_id']
        frame[f'dog_{name}'] = frame.gene_id.isin(ids)
    return frame


def dog_circular_shifts(d, n=10_000, p_col='min_p', count_col='n_snps', seed=20261013):
    """Original two-set circular shifts with the same RNG draw order.

    Batch shifted indices instead of rolling full label arrays per replicate.
    Sort indices to preserve the original order of scores in each mean.
    """
    covered = d.loc[d[count_col].gt(0)].copy()
    covered['snp_stratum'] = covered.groupby('chrom')[count_col].transform(
        lambda x: pd.qcut(x.rank(method='first'), min(5, len(x)), labels=False))
    covered['signal_rank'] = covered.groupby(['chrom', 'snp_stratum'])[p_col].transform(
        lambda x: x.rank(ascending=True, method='average', pct=True))
    groups = [x.sort_values('tss') for _, x in covered.groupby('chrom')]
    rng = np.random.default_rng(seed)
    out = {}
    for name, mask_col in [('all', 'dog'), ('significant', 'dog_significant')]:
        observed = float(covered.loc[covered[mask_col], 'signal_rank'].mean())
        offsets = np.array([[rng.integers(len(x)) for x in groups] for _ in range(n)])
        shifted_scores = []
        for j, group in enumerate(groups):
            positions = np.flatnonzero(group[mask_col].to_numpy())
            indices = np.sort((positions[None, :] + offsets[:, j, None]) % len(group), axis=1)
            shifted_scores.append(group.signal_rank.to_numpy()[indices])
        null = np.concatenate(shifted_scores, axis=1).mean(axis=1)
        out[name] = {'n': int(covered[mask_col].sum()),
                     'observed_mean_rank': observed, 'null_mean': float(null.mean()),
                     'empirical_one_sided_p': float((1 + (null <= observed).sum()) / (n + 1)),
                     'permutations': n}
    return covered, out
