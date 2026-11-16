import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.compare_psychiatric import read_psychiatric_gene_windows


class PsychiatricWindowFixtures(unittest.TestCase):
    def test_pgc_metadata_and_alternative_column_names(self):
        genes = pd.DataFrame({'gene_id': ['g1'], 'chrom': ['1'], 'tss': [500_001]})
        segments = {'1': [(0, 2_000_001, ('1', '1', 1, 0))]}
        for chrom_col in ['#CHROM', 'CHROM']:
            with self.subTest(chrom_col=chrom_col), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / 'pgc.tsv.gz'
                with gzip.open(path, 'wt') as handle:
                    handle.write('##genomeReference="GRCh37"\n##PVAL="P-value"\n')
                    handle.write(f'{chrom_col}\tPOS\tID\tPVAL\tBETA\n')
                    handle.write('1\t1000001\tright\t.02\t.1\n')
                    handle.write('1\t1\tleft\t.01\t.2\n')
                    handle.write('1\t20\tzero\t0\t.2\n')
                scores, audit = read_psychiatric_gene_windows(
                    path, genes, segments, chunksize=1,
                    columns={chrom_col: 'CHR', 'POS': 'BP', 'ID': 'SNP', 'PVAL': 'P'})
            self.assertEqual(scores.iloc[0].n_snps_500kb, 2)
            self.assertEqual(scores.iloc[0].min_p_500kb, .01)
            self.assertEqual(scores.iloc[0].lead_snp_500kb, 'left')
            self.assertEqual(audit['counts']['total'], 3)
            self.assertEqual(audit['counts']['p_zero'], 1)

    def test_unsorted_whitespace_liftover_boundaries_and_rejected_rows(self):
        # Forward/reverse blocks, an ambiguous block and a cross-chromosome map.
        segments = {'1': [(0, 2_000_001, ('1', '1', 1, 100)),
                          (2_000_001, 2_000_010, None),
                          (2_000_010, 2_000_020, ('1', '2', 1, -1_000_000))],
                    '2': [(0, 10, ('2', '2', -1, 1_000_009))]}
        genes = pd.DataFrame({'gene_id': ['g1', 'g2', 'empty'],
                              'chrom': ['1', '2', '3'],
                              'tss': [1_000_101, 1_000_010, 5_000_000]})
        rows = [
            'SNP  CHR BP P A1\n',
            'right 1 2000001 .2 A\n',  # inclusive +1 Mb endpoint after lift
            'narrow_right 1 1500001 .03 A\n',
            'reverse 2 1 .005 A\n',
            'left 1 1 .1 A\n',  # inclusive -1 Mb endpoint after lift
            'narrow_left 1 500001 .01 A\n',
            'cross 1 2000011 .006 A\n',
            'outside_narrow 1 500000 .001 A\n',
            'ambiguous 1 2000002 .00001 A\n',
            'unmapped 1 3000000 .00001 A\n',
            'sex X 1 .00001 A\n',
            'bad_position 1 1.5 .01 A\n',
            'zero 1 1 0 A\n',
            'bad_p 1 1 2 A\n',
            'missing_p 1 1 NA A\n',
            'infinite_p 1 1 inf A\n',
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'gwas.gz'
            with gzip.open(path, 'wt') as handle:
                handle.writelines(rows)
            scores, audit = read_psychiatric_gene_windows(
                path, genes, segments, chunksize=2)
        first, second, empty = scores.iloc[0], scores.iloc[1], scores.iloc[2]
        self.assertEqual(first.n_snps, 5)
        self.assertEqual(first.n_snps_500kb, 2)
        self.assertEqual(first.min_p, .001)
        self.assertEqual(first.min_p_500kb, .01)
        self.assertEqual(first.lead_pos_grch38_500kb, 500101)
        self.assertEqual(second.n_snps, 2)
        self.assertEqual(second.min_p, .005)
        self.assertEqual(second.lead_pos_grch38, 1000010)
        self.assertEqual(empty.n_snps, 0)
        self.assertTrue(pd.isna(empty.min_p))
        self.assertEqual(audit['counts']['total'], 15)
        self.assertEqual(audit['counts']['nonautosomal'], 1)
        self.assertEqual(audit['counts']['invalid_fields'], 5)
        self.assertEqual(audit['counts']['p_zero'], 1)
        self.assertEqual(audit['counts']['mapped'], 7)
        self.assertEqual(audit['counts']['ambiguous'], 1)
        self.assertEqual(audit['counts']['unmapped'], 1)
        self.assertEqual(audit['counts']['cross_chromosome_mapped'], 1)

    def test_rejects_file_without_valid_mapped_variants(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'empty.gz'
            with gzip.open(path, 'wt') as handle:
                handle.write('CHR SNP BP P\n1 bad 10 0\n')
            with self.assertRaisesRegex(ValueError, 'No valid mapped'):
                read_psychiatric_gene_windows(path, pd.DataFrame(), {})


if __name__ == '__main__':
    unittest.main()
