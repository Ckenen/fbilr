#!/usr/bin/env python
import os
import unittest
import filecmp

class TestFindBarcode(unittest.TestCase):
    def test_find_barcode(self):
        cmd = "./src/scripts/find_barcode.py -m tests/data/test.metrics.txt -s tests/data/test.stat.txt ./tests/data/barcodes.fa ./tests/data/test.fastq.gz"
        self.assertEqual(0, os.system(cmd))
        self.assertTrue(filecmp.cmp("tests/data/ref.metrics.txt", "tests/data/test.metrics.txt"))
        self.assertTrue(filecmp.cmp("tests/data/test.stat.txt", "tests/data/ref.stat.txt"))
    
if __name__ == "__main__":
    unittest.main()