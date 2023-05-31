#!/usr/bin/env python
import os
import time
import multiprocessing
import logging
from pygz import PigzFile
from fbilr import utils

DEBUG = True

CONFIDENCE_EDIT_DISTANCE = 3


class BarcodeFinder(object):
    def __init__(self):
        self.f_barcode = None
        self.f_reads = None
        self.threads = 1
        self.max_edit_distance = 5
        self.width = 200
        self.reads_per_batch = 10000 # 10k
        self.max_submitted_batch = self.threads * 2
        self.counter_of_summary = dict()
        self.f_matrix = None
        self.h_matrix = None
        self.f_summary = None
        self.splitted_fastq_dir = None
        self.h_fastqs = None
        self.ignore_read_name = False
    
    @staticmethod
    def extract_edge_sequences(reads, width):
        reference_sequences = []
        for read in reads:
            seq = read[1]
            length = len(seq)
            if length > width:
                seq_head = seq[:width]
                seq_tail = seq[-width:]
            else:
                seq_head = seq
                seq_tail = None
            reference_sequences.append([seq_head, seq_tail, length])
        return reference_sequences

    @staticmethod
    def _worker(references, barcodes):
        rows = []
        for seq_head, seq_tail, length in references:
            a = None  # alignment result: [loc1, loc2, ed]
            bc = None  # best barcode
            n_direction = 0  # 0: forward, 1: reverse
            n_location = 0  # 0: head 1: tail
            for bc_name, bc_seq_f, bc_seq_r in barcodes:
                tmp = utils.align(bc_seq_f, seq_head)
                if a is None or tmp[2] < a[2]:
                    a = tmp
                    bc, n_direction, n_location = bc_name, 0, 0
                    if a[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
                tmp = utils.align(bc_seq_r, seq_head)
                if tmp[2] < a[2]:
                    a = tmp
                    bc, n_direction, n_location = bc_name, 1, 0
                    if a[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
                if seq_tail is None:
                    continue
                tmp = utils.align(bc_seq_f, seq_tail)
                if tmp[2] < a[2]:
                    a = tmp
                    bc, n_direction, n_location = bc_name, 0, 1
                    if a[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break

                tmp = utils.align(bc_seq_r, seq_tail)
                if tmp[2] < a[2]:
                    a = tmp
                    bc, n_direction, n_location = bc_name, 1, 1
                    if a[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
            assert a
            direction = "F" if n_direction == 0 else "R"
            start, end = a[0], a[1]
            if n_location == 1:
                offset = length - len(seq_tail)
                start, end = start + offset, end + offset
            middle_position = int(length / 2)
            if end <= middle_position:
                location = "H"
            elif start >= middle_position:
                location = "T"
            else:
                location = "M"
            row = [bc, direction, location, start, end, a[2]]
            rows.append(row)
        return rows
        
    def prepare(self):
        self.max_submitted_batch = self.threads * 2
        # output for metrics
        if self.f_matrix.endswith(".gz"):
            self.h_matrix = PigzFile(self.f_matrix, "wt")
        else:
            self.h_matrix = open(self.f_matrix, "w+")
        # load barcodes
        self.barcodes = utils.load_barcodes(self.f_barcode)
        assert len(self.barcodes) > 0
        # summary
        self.counter_of_summary = dict()
        for barcode in self.barcodes:
            self.counter_of_summary[barcode[0]] = 0
        self.counter_of_summary["unclassified"] = 0
        # output fastq
        if self.splitted_fastq_dir:
            if not os.path.exists(self.splitted_fastq_dir):
                os.mkdir(self.splitted_fastq_dir)
            self.h_fastqs = dict()
            for barcode in self.barcodes:
                bc_name = barcode[0]
                self.h_fastqs[bc_name] = open(os.path.join(self.splitted_fastq_dir, bc_name + ".fastq"), "w+")
            self.h_fastqs["unclassified"] = open(os.path.join(self.splitted_fastq_dir, "unclassified.fastq"), "w+")
                        
    def finished(self):
        if self.h_matrix:
            self.h_matrix.close()
            
        if self.f_summary:
            with open(self.f_summary, "w+") as h:
                for k, v in self.counter_of_summary.items():
                    h.write("%s\t%d\n" % (k, v))
                    
        if self.h_fastqs:
            for h in self.h_fastqs.values():
                h.close()
    
    def run_pipeline(self):
        t1 = time.time()
        self.prepare()
        self.execute()
        self.finished()
        t2 = time.time()
        t = t2 - t1
        h = int(t / 3600)
        t = t % 3600
        m = t / 60
        s = t % 60
        logging.info("spent %dh%dm%.2fs" % (h, m, s))
        
    def execute(self):
        pool = multiprocessing.Pool(self.threads, maxtasksperchild=5)
        submitted_batch_list = []
        batch_loader = utils.load_batch(self.f_reads, self.reads_per_batch)
        for batch_id, reads in enumerate(batch_loader):    
            if DEBUG and batch_id >= 16:
                break
            reference_sequences = BarcodeFinder.extract_edge_sequences(reads, self.width)
            while True:
                # Priority process finished batch
                if len(submitted_batch_list) > 0:
                    item = submitted_batch_list[0]
                    if item[2].ready():
                        self.process_results(item[1], item[2].get())
                        submitted_batch_list.pop(0)
                        logging.info("Processed batch %d" % item[0])
                # Do not submit too many tasks
                if len(submitted_batch_list) < self.max_submitted_batch:
                    args = (reference_sequences, self.barcodes)
                    r = pool.apply_async(BarcodeFinder._worker, args)
                    submitted_batch_list.append([batch_id, reads, r])
                    logging.info("Submitted batch %d" % batch_id)
                    break
                else:
                    time.sleep(1)
        while len(submitted_batch_list) > 0:
            item = submitted_batch_list[0]
            if item[2].ready():
                self.process_results(item[1], item[2].get())
                submitted_batch_list.pop(0)
                logging.info("Processed batch %d" % item[0])
            else:
                time.sleep(1)
        pool.close()
        pool.join()
            
    def process_results(self, reads, rows):  
        assert len(reads) == len(rows)
        for read, row in zip(reads, rows):
            name, sequence, quality = read
            length = len(sequence)
            bc, ed = row[0], row[5]
            h_fastq = None                
            if ed <= self.max_edit_distance:
                self.counter_of_summary[bc] += 1
                if self.h_fastqs:
                    h_fastq = self.h_fastqs[bc]
            else:
                self.counter_of_summary["unclassified"] += 1
                if self.h_fastqs:
                    h_fastq = self.h_fastqs["unclassified"]
            if h_fastq:
                h_fastq.write("%s\n%s\n+\n%s\n" % (name, sequence, quality))
            if self.ignore_read_name:
                row1 = [".", length] # ignore read name
            else:
                row1 = [name.split()[0][1:], length]    
            row1.extend(row)
            self.h_matrix.write("\t".join(map(str, row1)) + "\n")
                        
            
class PairEndBarcodeFinder(BarcodeFinder):
    @staticmethod
    def _worker(references, barcodes):
        rows = []
        for seq_head, seq_tail, length in references:
            if length < len(seq_head) * 2:
                rows.append([".", ".", ".", -1, -1, -1, ".", ".", ".", -1, -1, -1])
                continue
            a_head = None # alignment result: [loc1, loc2, ed]
            a_tail = None
            bc_head = None # best barcode
            bc_tail = None
            n_direction_head = 0 # 0: forward, 1: reverse
            n_direction_tail = 0
            n_location_head = 0 # 0: head 1: tail
            n_location_tail = 0
            for bc_name, bc_seq_f, bc_seq_r in barcodes:
                tmp = utils.align(bc_seq_f, seq_head)
                if a_head is None or tmp[2] < a_head[2]:
                    a_head = tmp
                    bc_head, n_direction_head, n_location_head = bc_name, 0, 0
                    if a_head[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
                tmp = utils.align(bc_seq_r, seq_head)
                if tmp[2] < a_head[2]:
                    a_head = tmp
                    bc_head, n_direction_head, n_location_head = bc_name, 1, 0
                    if a_head[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
            for bc_name, bc_seq_f, bc_seq_r in barcodes:
                tmp = utils.align(bc_seq_f, seq_tail)
                if a_tail is None or tmp[2] < a_tail[2]:
                    a_tail = tmp
                    bc_tail, n_direction_tail, n_location_tail = bc_name, 0, 1
                    if a_tail[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
                tmp = utils.align(bc_seq_r, seq_tail)
                if tmp[2] < a_tail[2]:
                    a_tail = tmp
                    bc_tail, n_direction_tail, n_location_tail = bc_name, 1, 1
                    if a_tail[2] <= CONFIDENCE_EDIT_DISTANCE:
                        break
            assert a_head
            assert a_tail
            direction_head = "F" if n_direction_head == 0 else "R"
            direction_tail = "F" if n_direction_tail == 0 else "R"
            start_head, end_head = a_head[0], a_head[1]
            start_tail, end_tail = a_tail[0], a_tail[1]
            assert n_location_head == 0
            assert n_location_tail == 1
            offset = length - len(seq_tail)
            start_tail, end_tail = start_tail + offset, end_tail + offset
            row = [bc_head, direction_head, "H", start_head, end_head, a_head[2], 
                   bc_tail, direction_tail, "T", start_tail, end_tail, a_tail[2]]
            rows.append(row)
        return rows
    
    def prepare(self):
        self.max_submitted_batch = self.threads * 2
        # output for metrics
        if self.f_matrix.endswith(".gz"):
            self.h_matrix = PigzFile(self.f_matrix, "wt")
        else:
            self.h_matrix = open(self.f_matrix, "w+")
        # load barcodes
        self.barcodes = utils.load_barcodes(self.f_barcode)
        assert len(self.barcodes) > 0
        assert self.splitted_fastq_dir is None
        assert self.f_summary is None
    
    def execute(self):
        pool = multiprocessing.Pool(self.threads, maxtasksperchild=5)
        submitted_batch_list = []
        batch_loader = utils.load_batch(self.f_reads, self.reads_per_batch)
        for batch_id, reads in enumerate(batch_loader):    
            if DEBUG and batch_id >= 16:
                break
            reference_sequences = BarcodeFinder.extract_edge_sequences(reads, self.width)
            while True:
                # Priority process finished batch
                if len(submitted_batch_list) > 0:
                    item = submitted_batch_list[0]
                    if item[2].ready():
                        self.process_results(item[1], item[2].get())
                        submitted_batch_list.pop(0)
                        logging.info("Processed batch %d" % item[0])
                # Do not submit too many tasks
                if len(submitted_batch_list) < self.max_submitted_batch:
                    args = (reference_sequences, self.barcodes)
                    r = pool.apply_async(PairEndBarcodeFinder._worker, args)
                    submitted_batch_list.append([batch_id, reads, r])
                    logging.info("Submitted batch %d" % batch_id)
                    break
                else:
                    time.sleep(1)
        while len(submitted_batch_list) > 0:
            item = submitted_batch_list[0]
            if item[2].ready():
                self.process_results(item[1], item[2].get())
                submitted_batch_list.pop(0)
                logging.info("Processed batch %d" % item[0])
            else:
                time.sleep(1)
        pool.close()
        pool.join()
        
    def process_results(self, reads, rows):  
        assert len(reads) == len(rows)
        for read, row in zip(reads, rows):
            name, sequence, quality = read
            length = len(sequence)
            if self.ignore_read_name:
                row1 = [".", length] # ignore read name
            else:
                row1 = [name.split()[0][1:], length]    
            row1.extend(row)
            self.h_matrix.write("\t".join(map(str, row1)) + "\n")