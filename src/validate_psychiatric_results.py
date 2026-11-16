"""Independently check dog-gene minima/counts from the four psychiatric files.

Run after the Jupyter notebook: .venv/bin/python src/validate_psychiatric_results.py
This uses direct interval masks rather than the scoring reader's window search.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.compare_ocd import CHAIN, chain_segments, lift_positions
from src.compare_psychiatric import PSYCHIATRIC_GWAS, file_md5, read_psychiatric_chunks
from src.toga_gene_set import read_toga_gene_set

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=ROOT / 'notebook/results/ipynb/toga')
args = parser.parse_args()
output = args.output
_, selected = read_toga_gene_set(ROOT / 'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx',
    pd.read_csv(output / 'per_gene_ocd.csv'))
segments = chain_segments(CHAIN)
results = pd.read_csv(output / 'mann_whitney_results.csv')
expected_tests = (2 + len(PSYCHIATRIC_GWAS)) * 2 * 3
assert len(results) == expected_tests
primary = results.trait.eq('OCD')
assert primary.sum() == 6
np.testing.assert_allclose(results.loc[primary, 'holm_adjusted_p'],
    multipletests(results.loc[primary, 'one_sided_p'], method='holm')[1], rtol=1e-12, atol=1e-15)
assert results.loc[~primary, 'holm_adjusted_p'].isna().all()
ranked = pd.read_csv(output / 'mann_whitney_gene_ranks.csv')
validation = {}
for trait, spec in PSYCHIATRIC_GWAS.items():
    started = time.monotonic()
    frame = pd.read_csv(output / f'per_gene_{trait.lower()}.csv', dtype={'chrom': str})
    dog = frame.loc[frame.dog]
    assert set(dog.gene_id) == set(selected.human_gene_id)
    counts = np.zeros((len(dog), 2), dtype=np.int64)
    minima = np.full((len(dog), 2), np.inf)
    for chunk in read_psychiatric_chunks(
            ROOT / 'data/downloaded/gwas' / spec['filename'], columns=spec.get('columns')):
        for chrom, sub in chunk.groupby('CHR'):
            bp = pd.to_numeric(sub.BP, errors='coerce').to_numpy(dtype=float)
            p = pd.to_numeric(sub.P, errors='coerce').to_numpy(dtype=float)
            good = (np.isfinite(bp) & (bp >= 1) & (bp == np.floor(bp)) &
                    np.isfinite(p) & (p > 0) & (p <= 1))
            dest, pos, status = lift_positions(chrom, bp[good].astype(np.int64), segments)
            p = p[good]
            mapped = (status == 1) & (pos > 0)
            for i, gene in enumerate(dog.itertuples()):
                on_chrom = mapped & (dest == gene.chrom)
                if not on_chrom.any():
                    continue
                positions = pos[on_chrom]
                pvalues = p[on_chrom]
                for j, width in enumerate([1_000_000, 500_000]):
                    within = ((positions >= max(1, gene.tss - width)) &
                              (positions <= gene.tss + width))
                    counts[i, j] += int(within.sum())
                    if within.any():
                        minima[i, j] = min(minima[i, j], float(pvalues[within].min()))
    for j, (window, suffix) in enumerate([('1mb', ''), ('500kb', '_500kb')]):
        np.testing.assert_array_equal(counts[:, j], dog[f'n_snps{suffix}'])
        np.testing.assert_allclose(minima[:, j], dog[f'min_p{suffix}'], rtol=1e-12, atol=0)
        covered = frame.loc[frame[f'n_snps{suffix}'].gt(0)].copy()
        covered['expected_rank'] = covered[f'min_p{suffix}'].rank(method='average')
        stored = ranked.loc[ranked.trait.eq(trait) & ranked.window.eq(window)]
        check = covered.merge(stored[['gene_id', 'rank']], on='gene_id', validate='one_to_one')
        np.testing.assert_array_equal(check.expected_rank, check['rank'])
        for set_name, mask in [('significant', 'dog_significant'),
                               ('suggestive', 'dog_suggestive'), ('combined', 'dog')]:
            direct = mannwhitneyu(covered.loc[covered[mask], f'min_p{suffix}'],
                covered.loc[~covered.dog, f'min_p{suffix}'], alternative='less', method='asymptotic')
            saved = results.loc[results.trait.eq(trait) & results.window.eq(window) &
                                results.dog_set.eq(set_name)].iloc[0]
            np.testing.assert_allclose(direct.pvalue, saved.one_sided_p, rtol=1e-12, atol=1e-15)
    summary = json.loads((output / f'{trait.lower()}_summary.json').read_text())
    assert summary['scores_md5'] == file_md5(output / f'per_gene_{trait.lower()}.csv')
    # Reverse the sample point through all unique chain blocks.
    for sample in summary['gwas']['samples']:
        reverse = []
        for chrom, pieces in segments.items():
            for start, end, block in pieces:
                if block is None or block[1] != sample['grch38_chrom']:
                    continue
                _, _, sign, offset = block
                zero = sign * (sample['grch38'] - 1 - offset)
                if start <= zero < end:
                    reverse.append((chrom, int(zero + 1)))
        assert reverse == [(sample['chrom'], sample['grch37'])]
    validation[trait] = {'raw_dog_windows_checked': int(counts.size),
        'reverse_liftover_samples_checked': len(summary['gwas']['samples']),
        'genome_wide_ranks_and_six_mann_whitney_tests_verified': True,
        'seconds': time.monotonic() - started}
    print(f'{trait}: all {counts.size} dog windows, genome-wide ranks and six tests verified', flush=True)
validation['all_6_ocd_holm_adjustments_verified'] = True
validation['control_trait_holm_values_absent'] = True
(output / 'psychiatric_independent_validation.json').write_text(json.dumps(validation, indent=2))
