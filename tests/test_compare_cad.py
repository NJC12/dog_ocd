import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.compare_cad import read_cad_gene_windows


class CadWindowFixtures(unittest.TestCase):
    def test_grch38_windows_zero_p_and_invalid_rows(self):
        genes = pd.DataFrame({'gene_id': ['g1', 'g2'], 'chrom': ['1', '2'],
                              'tss': [1000000, 5000000]})
        lines = [
            'chromosome\tbase_pair_location\tp_value\trsid\teffect_allele\tother_allele\n',
            '1\t1\t0.1\trs1\tA\tG\n',
            '1\t500000\t0.2\trs2\tA\tG\n',
            '1\t1000000\t0.01\trs3\tA\tG\n',
            '1\t1500000\t0.03\tNA\tC\tT\n',
            '1\t1500001\t0\trs5\tA\tG\n',
            '1\t2000000\t0.5\trs6\tA\tG\n',
            '1\t2000001\t2\trs7\tA\tG\n',
            'X\t1\t0.01\trs8\tA\tG\n']
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'cad.tsv.gz'
            with gzip.open(path, 'wt') as handle:
                handle.writelines(lines)
            score, audit = read_cad_gene_windows(path, genes, chunksize=3)
        self.assertEqual(audit['total'], 8)
        self.assertEqual(audit['valid_autosomal'], 6)
        self.assertEqual(audit['invalid_fields'], 1)
        self.assertEqual(audit['nonautosomal'], 1)
        self.assertEqual(audit['p_zero'], 1)
        self.assertEqual(score.loc[0, 'n_snps'], 6)
        self.assertEqual(score.loc[0, 'min_p'], 0)
        self.assertEqual(score.loc[0, 'n_snps_500kb'], 3)
        self.assertEqual(score.loc[0, 'min_p_500kb'], .01)
        self.assertEqual(score.loc[1, 'n_snps'], 0)
        self.assertTrue(pd.isna(score.loc[1, 'min_p']))


if __name__ == '__main__':
    unittest.main()
