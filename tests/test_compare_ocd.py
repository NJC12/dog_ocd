import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.compare_ocd import chain_segments, lift_positions, read_genes, window_min


class ComparisonFixtures(unittest.TestCase):
    def test_tss_strand_and_protein_filter(self):
        text = ('chr1\tx\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "A";\n'
                'chr1\tx\tgene\t300\t400\t.\t-\t.\tgene_id "ENSG2.1"; gene_type "protein_coding"; gene_name "B";\n'
                'chr1\tx\tgene\t500\t600\t.\t+\t.\tgene_id "ENSG3.1"; gene_type "lncRNA"; gene_name "C";\n')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'genes.gtf.gz'
            with gzip.open(p,'wt') as f:f.write(text)
            g=read_genes(p)
        self.assertEqual(g.tss.tolist(),[100,400])

    def test_window_inclusive_overlap_and_empty(self):
        pos=np.array([1,100,200,300,301])
        p=np.array([.5,.4,.01,.3,.001])
        snp=np.array(['a','b','c','d','e'])
        self.assertEqual(window_min(200,pos,p,snp,100),(.01,'c',200,3))
        self.assertEqual(window_min(201,pos,p,snp,100),(.001,'e',301,3))
        self.assertEqual(window_min(999,pos,p,snp,100)[3],0)
        self.assertEqual(window_min(1,np.array([1]),np.array([.1]),np.array(['a']),100)[3],1)

    def test_500kb_window_is_nested_inside_1mb(self):
        pos=np.array([499999,500000,1000000,1500000,1500001,2000000])
        p=np.array([.001,.2,.3,.4,.0001,.00001])
        snp=np.array(list('abcdef'))
        wide=window_min(1000000,pos,p,snp,width=1000000)
        narrow=window_min(1000000,pos,p,snp,width=500000)
        self.assertEqual(wide, (.00001,'f',2000000,6))
        self.assertEqual(narrow, (.2,'b',500000,3))

    def test_liftover_one_based_and_ambiguity(self):
        chain=('chain 100 chr1 1000 + 100 110 chr1 1000 + 200 210 1\n10\n'
               '\nchain 50 chr1 1000 + 105 110 chr1 1000 + 300 305 2\n5\n'
               '\nchain 40 chr1 1000 + 300 310 chr1 1000 - 400 410 3\n10\n')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.chain';p.write_text(chain)
            seg=chain_segments(p)
        dest,mapped,status=lift_positions('1',np.array([100,101,105,106,110,111]),seg)
        self.assertEqual(status.tolist(),[0,1,1,2,2,0])
        self.assertEqual(mapped[:3].tolist(),[0,201,205])
        _,reverse_pos,reverse_status=lift_positions('1',np.array([301,310]),seg)
        self.assertEqual(reverse_status.tolist(),[1,1])
        self.assertEqual(reverse_pos.tolist(),[600,591])


if __name__=='__main__':unittest.main()
