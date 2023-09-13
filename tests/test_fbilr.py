#!/usr/bin/env python
import os
import unittest
import filecmp


class TestFBILR(unittest.TestCase):  
    
    def test_single_end_mode(self):
        cmd = "./src/scripts/fbilr -t 2 -b tests/data/20220708_GM12878.p5.fasta -o tests/data/20220708_GM12878.p5.matrix.tsv tests/data/20220708_GM12878.downsample.fastq.gz"
        self.assertEqual(0, os.system(cmd))
        self.assertTrue(filecmp.cmp("tests/data/20220708_GM12878.p5.matrix.tsv.ref", "tests/data/20220708_GM12878.p5.matrix.tsv"))
        
    def test_paired_end_mode(self):
        cmd = "./src/scripts/fbilr -t 2 -b tests/data/20220708_GM12878.p5_p7.fasta -o tests/data/20220708_GM12878.p5_p7.matrix.tsv tests/data/20220708_GM12878.downsample.fastq.gz"
        self.assertEqual(0, os.system(cmd))
        self.assertTrue(filecmp.cmp("tests/data/20220708_GM12878.p5_p7.matrix.tsv.ref", "tests/data/20220708_GM12878.p5_p7.matrix.tsv"))
    
    def test_multiple_barcode_list(self):
        cmd = "./src/scripts/fbilr -t 2 -b tests/data/20220708_GM12878.p5.fasta,tests/data/20220708_GM12878.p7.fasta -o tests/data/20220708_GM12878.p5.p7.matrix.tsv tests/data/20220708_GM12878.downsample.fastq.gz"
        self.assertEqual(0, os.system(cmd))
        self.assertTrue(filecmp.cmp("tests/data/20220708_GM12878.p5.p7.matrix.tsv.ref", "tests/data/20220708_GM12878.p5.p7.matrix.tsv"))
    
    
if __name__ == "__main__":
    unittest.main()