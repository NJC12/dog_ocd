import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from src.toga_gene_set import adjust_ocd_tests, apply_dog_labels, dog_circular_shifts, read_toga_gene_set


class TogaSelectionTests(unittest.TestCase):
    def test_holm_corrects_only_six_ocd_tests(self):
        results = pd.DataFrame([
            {'trait': trait, 'window': window, 'dog_set': group, 'one_sided_p': p}
            for trait in ['OCD', 'CAD']
            for window, values in [('1mb', [.001, .01, .02]), ('500kb', [.03, .04, .9])]
            for group, p in zip(['significant', 'suggestive', 'combined'], values)])
        results.loc[results.trait.eq('CAD'), 'one_sided_p'] = 1e-20
        adjusted = adjust_ocd_tests(results)
        np.testing.assert_allclose(adjusted.loc[adjusted.trait.eq('OCD'), 'holm_adjusted_p'],
                                   [.006, .05, .08, .09, .09, .9])
        self.assertTrue(adjusted.loc[adjusted.trait.eq('CAD'), 'holm_adjusted_p'].isna().all())
        self.assertTrue(adjusted.loc[adjusted.trait.eq('OCD'), 'holm_family_tests'].eq(6).all())
        with self.assertRaisesRegex(ValueError, 'six unique'):
            adjust_ocd_tests(results.iloc[1:])

    def test_original_workbook_restores_five_genes(self):
        root = Path(__file__).resolve().parents[1]
        original = pd.read_csv(root / 'notebook/results/per_gene_ocd.csv')
        dog, selected = read_toga_gene_set(root / 'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx', original)
        self.assertEqual(len(dog), 87)
        self.assertEqual(len(selected), 83)
        relabeled = apply_dog_labels(original, selected)
        self.assertEqual(set(relabeled.loc[relabeled.dog & ~original.dog, 'symbol']),
                         {'EMC2', 'PPA2', 'SH3GL3', 'PLXNA2', 'PROSER3'})
        self.assertEqual(int(relabeled.dog_significant.sum()), 43)
        self.assertEqual(int(relabeled.dog_suggestive.sum()), 40)
        columns = [c for c in original if not c.startswith('dog')]
        pd.testing.assert_frame_equal(original[columns], relabeled[columns])
        with self.assertRaisesRegex(ValueError, 'missing'):
            apply_dog_labels(original.loc[~original.symbol.eq('EMC2')], selected)

    def test_batched_shifts_match_original_roll_algorithm(self):
        frame = pd.DataFrame({'gene_id': list('abcdefghijklmno'),
                              'chrom': ['1'] * 9 + ['2'] * 6,
                              'tss': [9, 2, 1, 3, 4, 5, 6, 7, 8, 6, 5, 4, 3, 2, 1],
                              'min_p': [.1, .2, .1, .5, .3] * 3,
                              'n_snps': [1, 2, 3, 4, 5] * 3,
                              'dog': [True, False, True, False, False] * 3,
                              'dog_significant': [True, False, False, False, False] * 3})
        covered, actual = dog_circular_shifts(frame, n=73, seed=19)
        groups = [x.sort_values('tss') for _, x in covered.groupby('chrom')]
        rng = np.random.default_rng(19)
        for name, mask_col in [('all', 'dog'), ('significant', 'dog_significant')]:
            observed = covered.loc[covered[mask_col], 'signal_rank'].mean()
            null = []
            for _ in range(73):
                values = []
                for group in groups:
                    labels = np.roll(group[mask_col].to_numpy(), rng.integers(len(group)))
                    values.extend(group.signal_rank.to_numpy()[labels])
                null.append(np.mean(values))
            self.assertEqual(actual[name]['null_mean'], np.mean(null))
            self.assertEqual(actual[name]['empirical_one_sided_p'],
                             (1 + (np.array(null) <= observed).sum()) / 74)


if __name__ == '__main__':
    unittest.main()
