#!/usr/bin/env python
import sys
import optparse
import multiprocessing
from collections import defaultdict
import edlib
from Bio import SeqIO
from pygz import PigzFile


def align(seq, ref):
    return edlib.align(seq, ref, task='path', mode='HW')


def load_barcodes(path):
    barcodes = []
    if path.endswith(".gz"):
        f = PigzFile(path, "rt")
    else:
        f = open(path)
    for record in SeqIO.parse(f, "fasta"):
        name = record.name
        seq1 = str(record.seq)
        seq2 = str(record.seq.reverse_complement())
        barcodes.append((name, seq1, seq2))
    f.close()
    return barcodes


def load_reads(path):
    if path == "-":
        f = sys.stdin
    elif path.endswith(".gz"):
        f = PigzFile(path, "rt")
    else:
        f = open(path)
    for i, line in enumerate(f):
        j = i % 4
        if j == 0:
            k = line.find(" ")
            if k != -1:
                name = line[1:k]
            else:
                name = line[1:-1]
        elif j == 1:
            yield name, line[:-1]
    f.close()


def load_batch(path, count_per_batch):
    reads = None
    for read in load_reads(path):
        if reads is None:
            reads = [read]
        else:
            reads.append(read)
        if len(reads) >= count_per_batch:
            yield reads
            reads = None
    if reads is not None:
        yield reads
        reads = None


def find_barcode(seq_of_head, seq_of_tail, read_length, barcodes):
    CUTOFF = 3
    a = None  # align
    bc = 0  # barcode
    n_orient = 0  # 0: forward, 1: reverse
    n_loc = 0  # 0: head 1: tail
    for bc_name, bc_seq_forward, bc_seq_reverse in barcodes:
        tmp = align(bc_seq_forward, seq_of_head)
        if a is None or tmp["editDistance"] < a["editDistance"]:
            a = tmp
            bc, n_orient, n_loc = bc_name, 0, 0
        if a["editDistance"] <= CUTOFF:
            break

        tmp = align(bc_seq_reverse, seq_of_head)
        if tmp["editDistance"] < a["editDistance"]:
            a = tmp
            bc, n_orient, n_loc = bc_name, 1, 0
        if a["editDistance"] <= CUTOFF:
            break

        if seq_of_tail is None:
            continue

        tmp = align(bc_seq_forward, seq_of_tail)
        if tmp["editDistance"] < a["editDistance"]:
            a = tmp
            bc, n_orient, n_loc = bc_name, 0, 1
        if a["editDistance"] <= CUTOFF:
            break

        tmp = align(bc_seq_reverse, seq_of_tail)
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


def worker(task):
    barcodes, array = task
    results = []
    for seq1, seq2, length in array:
        r = find_barcode(seq1, seq2, length, barcodes)
        results.append(r)
    return results


USAGE = """find_barcode.py [options] barcodes.fasta reads.fastq.gz"""


def main():
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-t", "--threads", dest="threads",
                      type="int", default=1, 
                      help="Number of threads used.")
    parser.add_option("-m", "--metrics", dest="metrics", default=None, 
                      help="path of metrics file. [stdout]")
    parser.add_option("-s", "--stats", dest="stats", default=None, 
                      help="path of stats file. [stderr]")
    parser.add_option("-e", "--edit-istance", dest="ed", default=5)
    parser.add_option("-w", "--width", dest="width", default=200, 
                      help="Find barcode sequences within WIDTH of the boundary.")
    options, args = parser.parse_args()

    bc_fasta, in_fastq = args
    threads = options.threads
    width = options.width
    metrics = options.metrics
    stats = options.stats

    barcodes = load_barcodes(bc_fasta)
    counter = defaultdict(int)
    counter_ed_threshold = options.ed
    count_per_thread = 20000
    count_per_batch = threads * count_per_thread
    
    if metrics is None:
        mf = sys.stdout
    elif metrics.endswith(".gz"):
        mf = PigzFile(metrics, "wt")
    else:
        mf = open(metrics, "w+")

    # multi-threads
    if threads > 1:
        for reads in load_batch(in_fastq, count_per_batch):
            pool = multiprocessing.Pool(threads)
            tasks = []
            for i in range(0, len(reads), count_per_thread):
                j = min(i + count_per_thread, len(reads))
                array = []
                for k in range(i, j):
                    read = reads[k]
                    read_length = len(read[1])
                    seq_of_head = read[1][:width]
                    seq_of_tail = None
                    if len(read[1]) > width:
                        seq_of_tail = read[1][-width:]
                    array.append([seq_of_head, seq_of_tail, read_length])
                tasks.append([barcodes, array])
            res = pool.map(worker, tasks)
            pool.close()
            pool.join()

            i = 0
            for item in res:
                for x in item:
                    read = reads[i]
                    read_name = read[0]
                    read_length = len(read[1])
                    bc, orient, loc, loc1, loc2, ed = x
                    if ed <= counter_ed_threshold:
                        counter[bc] += 1
                    else:
                        counter["unclassified"] += 1
                    s = "\t".join(map(str, [read_name, read_length, bc, 
                                            orient, loc, loc1, loc2, ed]))
                    mf.write(s)
                    mf.write("\n")
                    i += 1
    # single thread
    else:
        for read in load_reads(in_fastq):
            read_name = read[0]
            read_length = len(read[1])
            seq_of_head = read[1][:width]
            seq_of_tail = None
            if len(read[1]) > width:
                seq_of_tail = read[1][-width:]
            bc, orient, loc, loc1, loc2, ed = find_barcode(
                seq_of_head, seq_of_tail, read_length, barcodes)
            if ed <= counter_ed_threshold:
                counter[bc] += 1
            else:
                counter["unclassified"] += 1
            s = "\t".join(map(str, [read_name, read_length, bc, 
                                    orient, loc, loc1, loc2, ed]))
            mf.write(s)
            mf.write("\n")

    mf.close()

    if stats is None:
        sf = sys.stderr
    else:
        sf = open(stats, "w+")
    for k in sorted(counter.keys()):
        v = counter[k]
        sf.write("%s\t%d\n" % (k, v))
    sf.close()

    exit(0)


if __name__ == '__main__':
    main()
