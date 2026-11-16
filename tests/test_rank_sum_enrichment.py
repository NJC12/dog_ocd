import unittest

import pandas as pd

from src.rank_sum_enrichment import rank_and_test


class RankSumFixtures(unittest.TestCase):
    def test_strongest_signal_ranking_combines_dog_categories_and_skips_uncovered(self):
        frame = pd.DataFrame({
            'gene_id': ['a', 'b', 'c', 'd', 'e', 'uncovered'],
            'symbol': ['A', 'B', 'C', 'D', 'E', 'F'],
            'chrom': ['1'] * 6,
            'dog': [True, True, False, False, False, False],
            'dog_significant': [True, False, False, False, False, False],
            'dog_suggestive': [False, True, False, False, False, False],
            'min_p': [.001, .01, .01, .2, .5, None],
            'n_snps': [2, 2, 2, 2, 2, 0],
        })
        ranks, results = rank_and_test(frame, 'OCD', '1mb', '', 1_000_000)
        by_set = {result['dog_set']: result for result in results}
        by_gene = ranks.set_index('gene_id')
        self.assertEqual(len(ranks), 5)
        self.assertEqual(by_gene.loc['a', 'rank'], 1)
        self.assertEqual(by_gene.loc['b', 'rank'], 2.5)
        self.assertEqual(by_gene.loc['c', 'rank'], 2.5)
        self.assertEqual(set(by_set), {'significant', 'suggestive', 'combined'})
        self.assertEqual([by_set[x]['n_dog'] for x in by_set], [1, 1, 2])
        self.assertTrue(all(result['n_other'] == 3 for result in results))
        self.assertAlmostEqual(by_set['significant']['u_statistic'], 3)
        self.assertAlmostEqual(by_set['suggestive']['u_statistic'], 2.5)
        self.assertAlmostEqual(by_set['combined']['u_statistic'], 5.5)
        self.assertAlmostEqual(by_set['combined']['probability_dog_stronger'], 5.5 / 6)
        self.assertAlmostEqual(by_set['combined']['one_sided_p'], 0.11806848511082646)

    def test_rejects_mismatched_dog_group(self):
        frame = pd.DataFrame({
            'gene_id': ['a', 'b'], 'symbol': ['A', 'B'], 'chrom': ['1', '1'],
            'dog': [True, False], 'dog_significant': [False, False],
            'dog_suggestive': [False, False], 'min_p': [.01, .1],
            'n_snps': [1, 1],
        })
        with self.assertRaisesRegex(ValueError, 'dog labels'):
            rank_and_test(frame, 'OCD', '1mb', '', 1_000_000)


if __name__ == '__main__':
    unittest.main()
