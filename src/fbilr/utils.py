#!/usr/bin/env python
import gzip
from Bio import SeqIO
import edlib
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
        barcodes.append((bc_name, bc_seq_f, bc_seq_r))
    f.close()
    return barcodes
    
def load_reads(path):
    if path.endswith(".gz"):
        f = PigzFile(path, "rt")
    else:
        f = open(path)
    read_name = None
    read_sequence = None
    read_quality = None
    for i, line in enumerate(f):
        j = i % 4
        if j == 0:
            read_name = line[:-1]
        elif j == 1:
            read_sequence = line[:-1]
        elif j == 3:
            read_quality = line[:-1]
            yield read_name, read_sequence, read_quality
    f.close()
    
def load_batch(path, reads_per_batch=100):
    reads = None
    for read in load_reads(path):
        if reads is None:
            reads = [read]
        else:
            reads.append(read)
        if len(reads) >= reads_per_batch:
            yield reads
            reads = None
    if reads is not None:
        yield reads
        
def align(query, reference):
    # len(query) >= len(reference)
    a = edlib.align(query, reference, task="locations", mode="HW")
    ed = a["editDistance"]
    start, end = a["locations"][0]
    end += 1
    return start, end, ed