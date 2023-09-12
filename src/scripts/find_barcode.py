#!/usr/bin/env python
import optparse
import logging
from fbilr.finders import BarcodeFinder
from fbilr.finders import PairEndBarcodeFinder



def main():
    
    logging.basicConfig(level=logging.INFO, format = '%(asctime)s - %(levelname)s: %(message)s')

    
    # Parameters
    
    parser = optparse.OptionParser(usage="find_barcode.py [options] barcodes.fa reads.fq.gz")
    
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1, 
                    help="Number of threads used. [%default]")
    parser.add_option("-e", "--max-edit-istance", dest="max_edit_distance", default=5,
                    help="Max edit distance of pass barcode. [%default]")
    parser.add_option("-c", "--confident-edit-distance", dest="confident_ed", default=0, 
                      help="Edit distance threshold for confident barcode. When edit distance <= threshold, stop to find other barcodes. [%default]")
    parser.add_option("-w", "--width", dest="width", default=200,
                    help="Find barcode sequences within WIDTH of the boundary. [%default]")
    parser.add_option("-m", "--matrix", dest="matrix", default=None,
                    help="Path of metrics file. [%default]")
    parser.add_option("-i", "--ignore-read-name", dest="ignore_read_name", action="store_true", default=False,
                      help="Ignore read name to save storage. [%default]")
    
    group = optparse.OptionGroup(parser, "For single-end barcode")
    group.add_option("-s", "--summary", dest="summary", default=None,
                    help="Path of stats file. [%default]")
    group.add_option("-d", "--outdir", dest="outdir",
                    help="Output directory for splitted fastq files. [%default]")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, "For paired-end barcode")
    group.add_option("-p", "--pair-end", dest="is_pair_end", action="store_true", default=False,
                    help="Find barcode at head and tail seperately. [%default]")
    parser.add_option_group(group)
    
    options, args = parser.parse_args()
    
    f_barcode = args[0]
    f_reads = args[1]
    f_matrix = options.matrix
    f_summary = options.summary
    outdir = options.outdir
    threads = options.threads
    max_edit_distance = options.max_edit_distance
    high_confident_ed = options.confident_ed
    width = options.width
    is_pair_end = options.is_pair_end
    ignore_read_name = options.ignore_read_name
    
    # Run
    
    if is_pair_end:
        obj = PairEndBarcodeFinder()
    else:
        obj = BarcodeFinder()
    obj.f_barcode = f_barcode
    obj.f_reads = f_reads
    obj.threads = threads
    obj.f_matrix = f_matrix
    obj.f_summary = f_summary
    obj.splitted_fastq_dir = outdir
    obj.max_edit_distance = max_edit_distance
    obj.high_confident_ed = high_confident_ed
    obj.width = width
    obj.ignore_read_name = ignore_read_name
    
    obj.run_pipeline()
    

if __name__ == "__main__":
    main()