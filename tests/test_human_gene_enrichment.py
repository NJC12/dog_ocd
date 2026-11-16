import unittest

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from src.compare_ocd import rank_and_permute
from src.human_gene_enrichment import circular_shift_gene_set, human_mann_whitney


class HumanGeneEnrichmentTests(unittest.TestCase):
    def test_mann_whitney_uses_complement_of_human_set_and_average_ties(self):
        frame = pd.DataFrame({
            'trait': ['OCD'] * 6, 'window': ['500kb'] * 6,
            'width_each_side_bp': [500_000] * 6,
            'gene_id': list('abcdef'), 'min_p': [.01, .01, .1, .2, .3, .4],
            'rank': [1.5, 1.5, 3, 4, 5, 6],
            'dog': [False, False, True, True, False, False],
        })
        result = human_mann_whitney(frame, {'a', 'b'}).iloc[0]
        expected = mannwhitneyu([.01, .01], [.1, .2, .3, .4],
                                alternative='less', method='asymptotic')
        self.assertEqual(result.n_human, 2)
        self.assertEqual(result.n_other, 4)  # Dog genes belong to the complement.
        self.assertEqual(result.u_statistic, 8)
        self.assertEqual(result.probability_human_stronger, 1)
        self.assertEqual(result.median_human_rank, 1.5)
        self.assertAlmostEqual(result.one_sided_p, expected.pvalue)
        self.assertEqual(result.holm_family_tests, 1)

    def test_circular_shifts_match_original_label_roll_with_same_seed(self):
        rng = np.random.default_rng(18)
        frame = pd.DataFrame({
            'gene_id': [f'g{i}' for i in range(72)],
            'chrom': np.repeat(['1', '2', '3'], 24),
            'tss': np.tile(np.arange(24), 3),
            'n_snps': rng.integers(1, 100, 72),
            'min_p': rng.uniform(.001, 1, 72),
        })
        frame['dog'] = frame.gene_id.isin({'g0', 'g2', 'g11', 'g28', 'g42'})
        frame['dog_significant'] = frame.dog
        old_covered, old_tests = rank_and_permute(frame, n=200)
        covered, stats = circular_shift_gene_set(frame,
            set(frame.loc[frame.dog, 'gene_id']), n=200)
        np.testing.assert_array_equal(covered.signal_rank, old_covered.signal_rank)
        np.testing.assert_array_equal(covered.snp_stratum, old_covered.snp_stratum)
        for name in ('n', 'permutations', 'empirical_one_sided_p'):
            self.assertEqual(stats[name], old_tests['all'][name])
        for name in ('observed_mean_rank', 'null_mean'):
            self.assertAlmostEqual(stats[name], old_tests['all'][name], places=14)

    def test_uncovered_human_gene_is_rejected(self):
        frame = pd.DataFrame({'gene_id': ['a', 'b'], 'chrom': ['1', '1'],
                              'tss': [1, 2], 'n_snps': [0, 1], 'min_p': [np.nan, .1]})
        with self.assertRaisesRegex(ValueError, 'covered'):
            circular_shift_gene_set(frame, {'a'}, n=10)


if __name__ == '__main__':
    unittest.main()
