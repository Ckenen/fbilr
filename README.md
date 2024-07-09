# FBILR: Find Barcode In Long Reads

## Description

FBILR is designed to find the best-matched barcode in long reads and report detailed information, such as direction, location, and edit distance. Since the barcode is likely to be located at one of the ends of the read (head or tail), and the read length is longer than 1,000 bp, FBILR restricts the search range to within 200 nt (`-w` option) of both ends to reduce the amount of computation and save time. Besides, FBILR can run in parallel (`-t` option).

In FBILR, edit distance represents the difference between barcode sequence and reference sequence, including mismatch, insertion, and deletion of bases. The edit distance is calculated by `edlib`.

For each barcode, FBILR searches the best-matched hits of the forward barcode in the read head, the reverse barcode in the read head, the forward barcode in the read tail, and the reverse barcode in the read tail, respectively. In single-end mode (`-m` option), report the minimum edit distance hit around all barcodes. In paired-end mode, report the minimum edit distance hit in read head and tail around all barcodes respectively.

![Schema](docs/schematic.png)

Here, we show the schema of the barcode that exists in a 100 nt read:

 * In case 1, the barcode exists in the head of the read with 0 edit distance (fully matched). 
 * In case 2, the barcode exists in the middle of the read with 2 edit distance (2 mismatch). 
 * In case 3, the barcode exists in the tail of the read with 3 edit distance (1 mismatch and 2 deletion).

Finally, the bar1 is the best-matched barcode in this read.

## Installation

    # 1. setuptools
    python setup.py test
    python setup.py install

    # 2. PyPI (TODO)
    pip install fbilr

    # 3. Manual installation
    export PYTHONPATH="`pwd`/src:${PYTHONPATH}"
    export PATH="`pwd`/src/scripts:${PATH}"


## Usage

The usage of FBILR is shown below:

    # Single-end
    fbilr -t 8 -w 200 -o matrix.tsv -b barcodes.fa reads.fq.gz
    fbilr -t 8 -w 200 -b barcodes.fa reads.fq.gz | pigz -p 8 -c > matrix.tsv.gz

    # Paired-end
    fbilr -t 8 -w 200 -m PE -o matrix.tsv -b barcodes.fa reads.fq.gz

    # Multiple barcode list
    fbilr -t 8 -w 200 -o matrix.tsv -b barcodes1.fa,barcodes2.fa reads.fq.gz

    # Ignore read name in output
    fbilr -t 8 -w 200 -i -o matrix.tsv -b barcodes.fa reads.fq.gz

    # Include read sequence and quality in output
    fbilr -t 8 -w 200 -q -o matrix.tsv -b barcodes.fa reads.fq.gz

    # Find barcode and split
    fbilr -t 8 -w 200 -q -b barcodes.fa reads.fq.gz | your_custom_split_script.py


## Output

The FBILR will output tab-delimited results that consist of multiple columns (shown as follows). In the results, one row corresponds to one read in the input FASTQ file. Each read can find an optimal barcode, even though the edit distance is large.

    column 1: read name, if the '-i' option is set, the value is '.'
    column 2: read length
    column 3: barcode name
    column 4: barcode orientation (F or R)
    column 5: barcode location (H, M or T)
    column 6: start in read (0-base, included)
    column 7: end in read (0-base, not included)
    column 8: edit distance
    column ...

    # Example:
    1b2e274b-9da7-4a5f-b40f-e6c36249d825    215     Bar4    R       T       172     196     0
    ed320d59-77c6-41ba-895d-f4fdba5855f2    249     Bar2    F       H       29      53      0
    9aa445f6-63b9-44e5-9b9c-43feea216b7a    492     Bar3    F       H       36      60      0
    3087cbe0-7b00-40ff-837c-4cc59cf7e7ff    280     Bar4    R       T       239     263     0
    15c53c45-ff43-4374-8716-049495d113aa    345     Bar4    F       H       27      50      3
    21c0fe8d-1725-42ba-b490-eec2cd6f76b3    408     Bar2    F       H       27      51      0
    90af744f-1367-493d-84e2-ca2375413e2d    551     Bar8    F       H       47      71      0

Column 3 to column 8 represent 1 hit (6 columns). The number of columns is flexible and depends on the number of barcode lists and mode. The structure of columns is: information columns (2) + hit columns (6 * N) + fastq columns (4, optional)

The number of columns in single-end mode is 2 + 6. The number of columns in paired-end mode is 2 + 6 * 2. If the `-q` option is set, an additional 4 columns (name, sequence, "+", quality) is append to the tail.

For 2 barcode lists, the number of columns in single-end mode is 2 + 6 * 2. The number of columns in paired-end mode is 2 + 6 * 2 * 2.

Number of barcode list|Mode|Include fastq|Number of column
:---:|:---:|:---:|:---:
1|Single-end|N|2 + 6 = 8
1|Single-end|Y|2 + 6 + 4 = 12
1|Paired-end|N|2 + 6 * 2 = 14
1|Paired-end|Y|2 + 6 * 2 + 4 = 18
2|Single-end|N|2 + 6 * 2 = 14
2|Single-end|Y|2 + 6 * 2 + 4 = 18
2|Paired-end|N|2 + 6 * 4 = 26
2|Paired-end|Y|2 + 6 * 4 + 4 = 30
3|Single-end|N|2 + 6 * 3 = 20
3|Single-end|Y|2 + 6 * 3 + 4 = 24
3|Paired-end|N|2 + 6 * 6 = 38
3|Paired-end|Y|2 + 6 * 6 + 4 = 42

## Splitting

### Example

    1. Demultiplexing XXX datasets.
    2. Demultiplexing XXX datasets.

## Packaging and distribute PyPI

    python -m build
    python3 -m twine upload --repository pypi dist/*

## Change logs

### 2023-09-13 (v1.2.0)

1. Added test for FBILR.
