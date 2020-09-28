### THESE FUNCTIONS WILL BE MOVED TO plot.py ###

import sys
import random

import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
import pysam

def plot_sample_len():
    fq_path = sys.argv[1]
    data = []
    for record in SeqIO.parse(fq_path, "fastq"):
        data.append(len(record))

    fig, ax1 = plt.subplots(1, 1, tight_layout=True)
    ax1.hist(data, bins=np.arange(0, max(data)+1+1)-0.5, color = "seagreen")
    plt.savefig("sample.len.png")

# def plot_by_length():
#     fpath = sys.argv[1]

#     tps = {}
#     total = {}
#     for line in open(fpath):
#         idx, vlen, found, *rest = line.strip('\n').split(' ')
#         vlen = int(vlen)
#         found = bool(found)
#         tps[vlen] = tps[vlen] + found if vlen in tps else 1
#         total[vlen] = total[vlen] + 1 if vlen in total else 1
#         # if found:
#         #     tps.append(vlen)
#         # total.append(vlen)

#     fig, ax1 = plt.subplots(1, 1, tight_layout=True)
#     R = [int(tps[k]/total[k]) for k in sorted(total.keys())]
#     print(R)
#     ax1.scatter(sorted(total.keys()), R)
#     # ax1.bar(sorted(total.keys()), [np.log(total[k]) for k in sorted(total.keys())], color="grey")
#     # ax1.bar(sorted(tps.keys()), [np.log(tps[k]) for k in sorted(tps.keys())], color="green")

#     # ax1.set_ylim(-1, 101)
#     plt.show()
#     # plt.savefig("tps_by_len.png")

def corr_len_abund():
    fpath = sys.argv[1]

    print("[W] Hardcoded cutoff: len<=500")

    lens = []
    abunds = []
    for line in open(fpath):
        idx, abund, l = line.split(',')
        abund, l = int(abund), int(l)
        if l<=500:
            lens.append(l)
            abunds.append(abund)

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    plt.figure(figsize=(8, 8))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # the scatter plot:
    ax_scatter.scatter(lens, abunds, color = "seagreen", alpha=0.75)

    # now determine nice limits by hand:
    binwidth = 1
    #max_len = max(lens)
    #max_ab = max(abunds)
    #ax_scatter.set_xlim((-500, max_len+500))
    #ax_scatter.set_ylim((-500, max_ab+500))

    # bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(lens, log=True, bins=np.arange(max(lens)-min(lens)+1+1+1)-0.5, color = "seagreen")
    ax_histy.hist(abunds, log=True, bins=np.arange(max(abunds)-min(abunds)+1+1+1)-0.5, color = "seagreen", orientation="horizontal")

    #ax_histx.set_xlim(ax_scatter.get_xlim())
    #ax_histy.set_ylim(ax_scatter.get_ylim())

    ax_scatter.grid(b=True, which='major', linestyle='--')
    ax_scatter.set_axisbelow(True)
    ax_histx.grid(b=True, which='major', linestyle='--')
    ax_histx.set_axisbelow(True)
    ax_histy.grid(b=True, which='major', linestyle='--')
    ax_histy.set_axisbelow(True)

    ax_scatter.set_xlabel("Lengths")
    ax_scatter.set_ylabel("Abundances")
    ax_histx.set_ylabel("Counts (log)")
    ax_histy.set_xlabel("Counts (log)")
    plt.savefig("corr.png")
    # plt.show()

def histo_len():
    fpath = sys.argv[1]

    lens = []
    min_l = float("inf")
    max_l = 0
    for line in open(fpath):
        idx, abund, l = line.split(',')
        l = int(l)
        lens.append(l)

    fig, ax1 = plt.subplots(1, 1, tight_layout=True)
    ax1.hist(lens, bins = np.arange(0,501)-0.5)

    ax1.set_xlim(0,510)
    ax1.set_xlabel("Length")
    ax1.set_ylabel("Count")

    plt.savefig("histo.len.png")


def plot_nal_by_err():
    bampath1 = sys.argv[1]
    bampath2 = sys.argv[2]

    fq_path = sys.argv[3]
    tot_strings = 0
    for record in SeqIO.parse(fq_path, "fastq"):
        tot_strings += 1
    # tot_strings = 11204873

    data1 = []
    bamfile1 = pysam.AlignmentFile(bampath1, "rb")
    for al in bamfile1.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        data1.append(bad_bases)
    print(min(data1), max(data1))

    data2 = []
    bamfile2 = pysam.AlignmentFile(bampath2, "rb")
    for al in bamfile2.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        data2.append(bad_bases)
    print(min(data2), max(data2))

    nbins = 7
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey = True, tight_layout=True, figsize=(13, 7))

    ax1.axhline(tot_strings, 0, 100, linewidth=1, color='r', label="Total specific")
    ax1.hist(data1, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", alpha=0.75, cumulative = True, histtype="step", label="Cumulative count")
    ax1.hist(data1, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", histtype="stepfilled", label="Count")
    ax1.set_title("HG00733")
    ax1.set_ylabel("# Primary Alignments")

    ax2.axhline(tot_strings, 0, 100, linewidth=1, color='r', label="Total specific")
    ax2.hist(data2, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", alpha=0.75, cumulative = True, histtype="step", label="Cumulative count")
    ax2.hist(data2, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", histtype="stepfilled", label="Count")
    ax2.set_title("NA19240")

    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xticks()
    plt.xlabel("# Errors")

    ax2.legend(loc=4, fancybox=True, shadow=True)
    # plt.savefig("iden.2.png")
    plt.show()

if __name__ == "__main__":
    # plot_sample_len()
    # corr_len_abund()
    plot_nal_by_err()
    # plot_nal_by_id()
    # histo_len()
