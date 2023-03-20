#!/usr/bin/env python
import sys
import os
import optparse
import multiprocessing
import time
from collections import defaultdict
import edlib
from Bio import SeqIO
from pygz import PigzFile



USAGE = """

    find_barcode.py [options] barcodes.fasta reads.fastq.gz"""


class FindBarcode(object):
    def __init__(self):
        self.barcode_fasta = None
        self.read_fastq = None
        self.threads = 1
        self.width = 200
        self.metrics_file = None
        self.stats_file = None
        self.outdir = None
        self.max_edit_distance = 5
        self.reads_per_thread = 1000 # 100k

        self.barcodes = None
        self.fhs = None
        self.mf = None
        self.counter = defaultdict(int)
        
        self.init_parameters()
        self.execute()
        
    
    def init_parameters(self):
        parser = optparse.OptionParser(usage=USAGE)
        parser.add_option("-t", "--threads", dest="threads", type="int", default=1, 
                        help="Number of threads used.")
        parser.add_option("-e", "--edit-istance", dest="ed", default=5, 
                          help="Max edit distance of pass barcode.")
        parser.add_option("-w", "--width", dest="width", default=200, 
                        help="Find barcode sequences within WIDTH of the boundary.")
        parser.add_option("-m", "--metrics", dest="metrics", default=None, 
                        help="path of metrics file. [stdout]")
        parser.add_option("-s", "--stats", dest="stats", default=None, 
                        help="path of stats file. [stderr]")
        parser.add_option("-d", "--outdir", dest="outdir",
                        help="Output directory for spitted fastq.")
        options, args = parser.parse_args()
            
        self.barcode_fasta, self.read_fastq = args
        self.threads = options.threads
        self.width = options.width
        self.metrics_file = options.metrics
        self.stats_file = options.stats
        self.outdir = options.outdir
        self.max_edit_distance = options.ed

    def load_barcodes(self):
        barcodes = []
        if self.barcode_fasta.endswith(".gz"):
            f = PigzFile(self.barcode_fasta, "rt")
        else:
            f = open(self.barcode_fasta)
        for record in SeqIO.parse(f, "fasta"):
            bc_name = record.name
            bc_seq_f = str(record.seq)
            bc_seq_r = str(record.seq.reverse_complement())
            barcodes.append((bc_name, bc_seq_f, bc_seq_r))
        f.close()
        self.barcodes = barcodes

    def load_reads(self):
        if self.read_fastq == "-":
            f = sys.stdin
        elif self.read_fastq.endswith(".gz"):
            f = PigzFile(self.read_fastq, "rt")
        else:
            f = open(self.read_fastq)
        name = None
        sequence = None
        quality = None
        for i, line in enumerate(f):
            j = i % 4
            if j == 0:
                name = line[:-1]
            elif j == 1:
                sequence = line[:-1]
            elif j == 3:
                quality = line[:-1]
                yield name, sequence, quality
        f.close()
    
    def load_batch(self):
        reads = None
        for read in self.load_reads():
            if reads is None:
                reads = [read]
            else:
                reads.append(read)
            if len(reads) >= self.reads_per_thread:
                yield reads
                reads = None
        if reads is not None:
            yield reads
        
    def process_results(self, reads, rows):
        assert len(reads) == len(rows)
        for read, row in zip(reads, rows):
            name, sequence, quality = read
            bc, orient, loc, loc1, loc2, ed = row
            length = len(sequence)
            if ed <= self.max_edit_distance:
                s = "pass"
                self.counter[bc] += 1
            else:
                s = "fail"
                self.counter["unclassified"] += 1
            if self.fhs:
                h = self.fhs[(bc, "%s%s" % (orient, loc), s)]
                h.write("%s\n%s\n+\n%s\n" % (name, sequence, quality))
            row1 = [name.split()[0][1:], length]
            row1.extend(row)
            self.mf.write("\t".join(map(str, row1)) + "\n")

    def execute(self):
        self.load_barcodes()
        
        # output for splitted fastqs        
        if self.outdir:
            if not os.path.exists(self.outdir):
                os.mkdir(self.outdir)
            fhs = dict()
            for item in self.barcodes:
                bc_name = item[0]
                bc_dir = "%s/%s" % (self.outdir, bc_name)
                if not os.path.exists(bc_dir):
                    os.mkdir(bc_dir)
                for s in ["FH", "FT", "FM", "RH", "RT", "RM"]:
                    fhs[(bc_name, s, "pass")] = open("%s/%s_pass.fastq" % (bc_dir, s), "w+")
                    fhs[(bc_name, s, "fail")] = open("%s/%s_fail.fastq" % (bc_dir, s), "w+")
            self.fhs = fhs
    
        # output for metrics
        if self.metrics_file is None:
            self.mf = sys.stdout
        elif self.metrics_file.endswith(".gz"):
            self.mf = PigzFile(self.metrics_file, "wt")
        else:
            self.mf = open(self.metrics_file, "w+")
            
        # find barcode
        pool = multiprocessing.Pool(self.threads, maxtasksperchild=5)
        array = []
        for i, reads in enumerate(self.load_batch()):
            if i >= 8:
                break
            while True:
                if len(array) > 0:
                    if array[0][2].ready():
                        self.process_results(array[0][1], array[0][2].get())
                        array.pop(0)
                if len(array) < self.threads * 2:
                    args = (reads, self.barcodes, self.width)
                    r = pool.apply_async(FindBarcode._worker, args)
                    array.append([i, reads, r])
                    break
                else:
                    time.sleep(1)
        while len(array) > 0:
            if array[0][2].ready():
                self.process_results(array[0][1], array[0][2].get())
                array.pop(0)
            else:
                time.sleep(1)
        pool.close()
        pool.join()
        
        # close handles
        if self.fhs:
            for fh in fhs.values():
                fh.close()
        if self.mf:
            self.mf.close()
            
        # stats
        if self.stats_file is None:
            sf = sys.stderr
        else:
            sf = open(self.stats_file, "w+")
        for k in sorted(self.counter.keys()):
            v = self.counter[k]
            sf.write("%s\t%d\n" % (k, v))
        sf.close()
        
    @staticmethod
    def _align(seq, ref):
        return edlib.align(seq, ref, task='locations', mode='HW')

    @staticmethod
    def _find_barcode(seq_of_head, seq_of_tail, read_length, barcodes):
        CUTOFF = 3
        a = None  # align
        bc = 0  # barcode
        n_orient = 0  # 0: forward, 1: reverse
        n_loc = 0  # 0: head 1: tail
        for bc_name, bc_seq_f, bc_seq_r in barcodes:
            tmp = FindBarcode._align(bc_seq_f, seq_of_head)
            if a is None or tmp["editDistance"] < a["editDistance"]:
                a = tmp
                bc, n_orient, n_loc = bc_name, 0, 0
            if a["editDistance"] <= CUTOFF:
                break

            tmp = FindBarcode._align(bc_seq_r, seq_of_head)
            if tmp["editDistance"] < a["editDistance"]:
                a = tmp
                bc, n_orient, n_loc = bc_name, 1, 0
            if a["editDistance"] <= CUTOFF:
                break

            if seq_of_tail is None:
                continue

            tmp = FindBarcode._align(bc_seq_f, seq_of_tail)
            if tmp["editDistance"] < a["editDistance"]:
                a = tmp
                bc, n_orient, n_loc = bc_name, 0, 1
            if a["editDistance"] <= CUTOFF:
                break

            tmp = FindBarcode._align(bc_seq_r, seq_of_tail)
            if tmp["editDistance"] < a["editDistance"]:
                a = tmp
                bc, n_orient, n_loc = bc_name, 1, 1
            if a["editDistance"] <= CUTOFF:
                break
        assert a

        orient = "F" if n_orient == 0 else "R"
        loc1, loc2 = a["locations"][0]
        loc2 += 1
        if n_loc == 1:
            offset = read_length - len(seq_of_tail)
            loc1, loc2 = loc1 + offset, loc2 + offset
        md = int(read_length / 2)
        if loc2 <= md:
            loc = "H"
        elif loc1 >= md:
            loc = "T"
        else:
            loc = "M"
        ed = a["editDistance"]
        return bc, orient, loc, loc1, loc2, ed

    @staticmethod
    def _worker(reads, barcodes, width):  
        rows = []
        for read in reads:
            seq = read[1]
            length = len(seq)
            seq_head = seq
            seq_tail = None
            if length > width:
                seq_head = seq[:width]
                seq_tail = seq[-width:]
            row = FindBarcode._find_barcode(seq_head, seq_tail, length, barcodes)
            rows.append(row)
        return rows

    
    
if __name__ == '__main__':
    FindBarcode()
