#!/usr/bin/env python
import sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pygz import GzipFile


MAX_LENGTH = 5000
BIN_WIDTH = 100
BIN_COUNT = int(MAX_LENGTH / BIN_WIDTH)
MAX_EDIT_DISTANCE = 5


def plot_edit_distance_distribution(rows, outfile):
    counter = defaultdict(int)
    for row in rows:
        counter[row[6]] += 1
    ys = np.zeros(15)
    for i in range(len(ys)):
        ys[i] = counter[i]
    xs = np.arange(len(ys))
    ys = ys * 100 / len(rows)

    plt.figure(figsize=(4, 3))
    plt.bar(xs, ys, color="C0", lw=1, edgecolor="black")
    plt.xticks(xs)
    plt.xlim(min(xs) - 0.5, max(xs) + 0.5)
    plt.ylim(0, 100)
    plt.xlabel("Edit distance")
    plt.ylabel("Percentage (%)")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def plot_orient_loc_distribution(rows, outfile):
    counter = dict()
    for row in rows:
        d = row[6]
        if d not in counter:
            counter[d] = defaultdict(int)
        counter[d][(row[2], row[3])] += 1

    fig, axs = plt.subplots(4, 5, figsize=(12, 10), sharex=True, sharey=True)

    for i in range(20):
        m = int(i / 5)
        n = i % 5
        if i in counter:
            d = counter[i]
        else:
            d = defaultdict(int)
        keys = [("F", "H"), ("F", "T"), ("R", "H"), ("R", "T"), ]
        xticks = ["FH", "FT", "RH", "RT"]
        counts = []
        for key in keys:
            count = d[key]
            counts.append(count)
        ys = np.array(counts) * 100 / sum(counts)
        xs = np.arange(len(ys))
        plt.sca(axs[m][n])
        plt.title("ED=%d, N=%d" % (i, sum(counts)))
        plt.bar(xs, ys, color="C0", edgecolor="black", lw=1)
        plt.ylim(0, 100)
        plt.xticks(xs, xticks)
        plt.ylabel("Percentage (%)")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def plot_length_distribution(rows, outfile):
    array = []
    for i in range(MAX_EDIT_DISTANCE):
        array.append(np.zeros(BIN_COUNT))
    counts = np.zeros(BIN_COUNT)

    for row in rows:
        i = int(row[0] / BIN_WIDTH)
        if 0 <= i < len(counts):
            counts[i] += 1
            if row[6] < MAX_EDIT_DISTANCE:
                array[row[6]][i] += 1

    counts = counts / BIN_WIDTH
    for i in range(len(array)):
        array[i] = array[i] / BIN_WIDTH

    xs = (np.arange(BIN_COUNT) + 0.5) * BIN_WIDTH

    fig, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

    plt.sca(axs[0])
    plt.plot(xs, counts, label="All")
    plt.plot(xs, sum(array), label="ED<=%d" % MAX_EDIT_DISTANCE)
    plt.ylabel("Read numbers")
    plt.legend()
    
    ratios = sum(array) / counts
    plt.sca(axs[1])
    plt.plot(xs, ratios)
    plt.ylim(0, 1)
    plt.ylabel("Ratio of barcode/total")
    plt.xlabel("Read length (nt)")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def plot_length_edit_distance(rows, outfile):
    counter = dict()
    for row in rows:
        i = int(row[0] / BIN_WIDTH)
        if i not in counter:
            counter[i] = defaultdict(int)
        counter[i][row[6]] += 1
    xs = (np.arange(BIN_COUNT) + 0.5) * BIN_WIDTH

    values = np.zeros(BIN_COUNT)
    for i in range(BIN_COUNT):
        if i in counter:
            v = sum([x for x in counter[i].values()])
        else:
            v = 0
        values[i] = v

    plt.figure(figsize=(6, 4))
    bottoms = np.zeros(BIN_COUNT)
    for ed in range(10):
        ys = np.zeros(BIN_COUNT)
        for i in range(BIN_COUNT):
            if i in counter:
                y = counter[i][ed]
            else:
                y = 0
            ys[i] = y
        ys = ys / values
        plt.bar(xs, ys, width=BIN_WIDTH, bottom=bottoms, label=str(ed))
        bottoms += ys
    plt.bar(xs, 1 - bottoms, width=BIN_WIDTH,
            bottom=bottoms, color="lightgrey", label=">9")
    plt.xlabel("Read length (nt)")
    plt.ylabel("Proportion")
    plt.xlim(0, MAX_LENGTH)
    plt.ylim(0, 1)
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    

def plot_location(rows, outfile):
    counter_f_h = defaultdict(int)
    counter_f_t = defaultdict(int)
    counter_r_h = defaultdict(int)
    counter_r_t = defaultdict(int)

    for row in rows:
        if row[6] > MAX_EDIT_DISTANCE: # edit distance
            continue
        orient = row[2]
        length = row[0]
        start = row[4]
        end = row[5]
        dis1 = start
        dis2 = length - end
        if orient == "F":
            counter_f_h[dis1] += 1
            counter_f_t[dis2] += 1
        else:
            counter_r_h[dis1] += 1
            counter_r_t[dis2] += 1

    fig, axs = plt.subplots(2, 2, figsize=(10, 6))

    plt.sca(axs[0][0])
    plt.title("Barcodes in forward")
    xs = np.arange(0, 101)
    ys = [counter_f_h[x] for x in xs]
    t = sum(counter_f_h.values())
    ys = np.array(ys) * 100 / t
    plt.bar(xs, ys, width=1)
    plt.xlim(min(xs) - 1, max(xs) + 1)
    plt.ylabel("Percentage (%)")
    plt.xlabel("Distance to head end")
    plt.tight_layout()

    plt.sca(axs[0][1])
    plt.title("Barcodes in reverse")
    xs = np.arange(-100, 1)
    ys = [counter_f_t[-x] for x in xs]
    t = sum(counter_f_t.values())
    ys = np.array(ys) * 100 / t
    plt.bar(xs, ys, width=1)
    plt.xlim(min(xs) - 1, max(xs) + 1)
    plt.xlabel("Distance to tail end")
    plt.tight_layout()

    plt.sca(axs[1][0])
    plt.title("Barcodes in forward")
    xs = np.arange(0, 101)
    ys = [counter_r_h[x] for x in xs]
    t = sum(counter_r_h.values())
    ys = np.array(ys) * 100 / t
    plt.bar(xs, ys, width=1)
    plt.xlim(min(xs) - 1, max(xs) + 1)
    plt.ylabel("Percentage (%)")
    plt.xlabel("Distance to head end")
    plt.tight_layout()

    plt.sca(axs[1][1])
    plt.title("Barcodes in reverse")
    xs = np.arange(-100, 1)
    ys = [counter_r_t[-x] for x in xs]
    t = sum(counter_r_t.values())
    ys = np.array(ys) * 100 / t
    plt.bar(xs, ys, width=1)
    plt.xlim(min(xs) - 1, max(xs) + 1)
    plt.xlabel("Distance to tail end")
    plt.tight_layout()

    plt.savefig(outfile, dpi=300)


def main():
    infile, prefix = sys.argv[1:]

    rows = []
    if infile.endswith(".gz"):
        f = GzipFile(infile, "rt")
    else:
        f = open(infile)
    for line in f:
        row = line.strip("\n").split("\t")
        row[1] = int(row[1])
        row[5] = int(row[5])
        row[6] = int(row[6])
        row[7] = int(row[7])
        rows.append(row[1:]) # remove read name
    f.close()
            
    plot_edit_distance_distribution(rows, prefix + ".ed_dis.pdf")
    plot_orient_loc_distribution(rows, prefix + ".orient_loc.pdf")
    plot_length_distribution(rows, prefix + ".len.pdf")
    plot_length_edit_distance(rows, prefix + ".len_ed.pdf")
    plot_location(rows, prefix + ".loc_dis.pdf")


if __name__ == '__main__':
    main()
