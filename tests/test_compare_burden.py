import unittest

import pandas as pd

from src.compare_burden import analyze


class BurdenTests(unittest.TestCase):
    def test_id_join_missing_genes_ties_and_u_direction(self):
        genes = pd.DataFrame(dict(gene_id=list('abcde'), symbol=list('ABCDE'),
                                  chrom=['1']*5, dog=[True, True, False, False, True],
                                  dog_significant=[True, False, False, False, True],
                                  dog_suggestive=[False, True, False, False, False]))
        burden = pd.DataFrame({'Gene ID': list('dcba'), 'P SKAT-O': [.9, .2, .2, .01]})
        ranked, results = analyze(genes, burden)
        self.assertEqual(ranked.set_index('gene_id')['rank'].to_dict(), dict(a=1, b=2.5, c=2.5, d=4))
        combined = results.set_index('dog_set').loc['combined']
        self.assertEqual(combined.n_dog, 2)
        self.assertEqual(combined.n_other, 2)
        self.assertEqual(combined.u_statistic, 3.5)
        self.assertEqual(combined.probability_dog_stronger, .875)


if __name__ == '__main__':
    unittest.main()
