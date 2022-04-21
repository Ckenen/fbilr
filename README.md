# FBILR: Find Barcode In Long Reads

## Description

FBILR is designed to find the best-matched barcode on reads and report detailed information (orientation, location, and edit distance). Since the barcode is likely to be located at one of the ends of the read (head or tail), and the read length is longer than 1,000 bp, FBILR restricts the search range to within 200bp of both ends to reduce the amount of computation and save time. Besides, FBILR is able to run in parallel and the order of output is accordant to input.

## Usage

The usage of FBILR is shown below:

    find_barcodes.py -t 8 -w 200 -m metrics.txt -s stats.txt barcodes.fasta reads.fastq.gz
    # or
    find_barcodes.py -t 8 -w 200 -s stats.txt barcodes.fasta reads.fastq.gz | gzip -c > metrics.txt.gz
    # or 
    pigz -p 4 -d -c reads.fastq.gz | find_barcodes.py -t 8 -w 200 -s stats.txt barcodes.fasta - | pigz -p 4 -c > metrics.txt.gz

After that, you can visualize the result with the following command:

    plot_barcode_detail.py -m metrics.txt -p outdir/out
    # or 
    pigz -p 4 -d -c metrics.txt.gz | plot_barcode_detail.py -p outdir/out

## Output files

The find_barcode.py script will output 2 files (`metrics.txt` and `stats.txt`).

The `metrics.txt` is a tab-delimited file that consists of 8 columns (shown as follows). In this file, one row corresponds to one read in the input FASTQ file. Each read can find an optimal barcode, even though the edit distance is large (edit distance represents the difference between barcode sequence and reference sequence, including mismatch, insertion, and deletion of bases).

    column 1: read name
    column 2: read length
    column 3: barcode name
    column 4: barcode orientation (F or R)
    column 5: barcode location (H, M or T)
    column 6: start in read (0-base, included)
    column 7: end in read (0-base, not included)
    column 8: edit distance

    # Example:
    1b2e274b-9da7-4a5f-b40f-e6c36249d825    215     Bar4    R       T       172     196     0
    ed320d59-77c6-41ba-895d-f4fdba5855f2    249     Bar2    F       H       29      53      0
    9aa445f6-63b9-44e5-9b9c-43feea216b7a    492     Bar3    F       H       36      60      0
    3087cbe0-7b00-40ff-837c-4cc59cf7e7ff    280     Bar4    R       T       239     263     0
    15c53c45-ff43-4374-8716-049495d113aa    345     Bar4    F       H       27      50      3
    21c0fe8d-1725-42ba-b490-eec2cd6f76b3    408     Bar2    F       H       27      51      0
    90af744f-1367-493d-84e2-ca2375413e2d    551     Bar8    F       H       47      71      0

The `stats.txt` is a tab-delimited file that consists of 2 columns (shown ad follows). Whether a barcode exists on a read is determined by edit distance (similarity). If the edit distance is small enough (-e option), we can confidently judge that the barcode exists, otherwise does not exist (unclassified). This file statistics the count of each barcode.

    column 1: barcode name
    column 2: count

    # Example: 
    Bar1	783154
    Bar2	236937
    Bar3	1579564
    Bar4	1266932
    Bar40	1983196
    Bar5	2876236
    Bar6	1571845
    Bar7	1663693
    Bar8	5781377
    Bar9	720325
    unclassified	2578447


## Schema

Here, we show the schema of the barcode that exists in the read (100 nt):

![Schema](src/schema.png)

In case 1, the barcode exists in the head of the read with 0 edit distance (fully matched). 

In case 2, the barcode exists in the middle of the read with 2 edit distance (2 mismatch). 

In case 3, the barcode exists in the tail of the read with 3 edit distance (1 mismatch and 2 deletion).

Finally, the bar1 is the best-matched barcode in this read.
