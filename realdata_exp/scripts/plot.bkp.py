### THESE FUNCTIONS - IF NEEDED - MUST BE MOVED TO plot.py ###

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



if __name__ == "__main__":
    # plot_sample_len()
    histo_len()
