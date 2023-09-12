#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO
import edlib
import pysam
from pygz import PigzFile


def load_barcodes(path):
    barcodes = []
    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path)
    for record in SeqIO.parse(f, "fasta"):
        bc_name = record.name
        bc_seq_f = str(record.seq)
        bc_seq_r = str(record.seq.reverse_complement())
        barcodes.append([bc_name, bc_seq_f, bc_seq_r])
    f.close()
    return barcodes
    
    
def load_reads(path):
    count = 0 # loaded read number
    if path == "-":
        f = sys.stdin
    if path.endswith(".gz"):
        f = PigzFile(path, "rt")
    else:
        f = open(path)
    name = None
    seq = None
    qua = None
    for i, line in enumerate(f):
        j = i % 4
        if j == 0:
            name = line[:-1]
        elif j == 1:
            seq = line[:-1]
        elif j == 3:
            qua = line[:-1]
            yield [name, seq, qua]
            count += 1
    f.close()
    
    
def load_batch(path, reads_per_batch=100):
    reads = list()
    with pysam.FastqFile(path) as f:
        for read in f:
            reads.append([read.name, read.sequence, read.quality])
            if len(reads) >= reads_per_batch:
                yield reads
                reads = list()
    if len(reads) > 0:
        yield reads
    
    # reads = None
    # for read in load_reads(path):
    #     if reads is None:
    #         reads = [read]
    #     else:
    #         reads.append(read)
    #     if len(reads) >= reads_per_batch:
    #         yield reads
    #         reads = None
    # if reads is not None:
    #     yield reads


def cut_edge_sequence(seq, width):
    head, tail = "", ""
    if len(seq) > 0 and width > 0:
        width = min(width, len(seq))
        head, tail = seq[:width], seq[-width:]
    return [len(seq), head, tail]

        
def align(query, reference):
    # assert len(query) <= len(reference)
    a = edlib.align(query, reference, task="locations", mode="HW")
    ed = a["editDistance"]
    start, end = a["locations"][0]
    end += 1
    return start, end, ed