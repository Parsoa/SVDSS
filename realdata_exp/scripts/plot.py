import sys

from Bio import SeqIO
import pysam

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_samplelengths():
    fq_path1 = sys.argv[1]
    fq_path2 = sys.argv[2]
    out_path = sys.argv[3]

    data1 = []
    for record in SeqIO.parse(fq_path1, "fastq"):
        l = len(record)
        data1.append(l)

    data2 = []
    for record in SeqIO.parse(fq_path2, "fastq"):
        l = len(record)
        data2.append(l)

    data = {"Sample" : ["HG00733"]*len(data1) + ["NA19240"]*len(data2),
            "Lengths" : data1 + data2}
    df = pd.DataFrame(data = data)

    dplot = sns.displot(df, x="Lengths", hue="Sample", fill=True, legend=False)
    dplot.ax.legend(labels=["HG00733", "NA19240"], loc=2)

    plt.savefig(out_path)
    # plt.show()

def plot_correlation():
    fq_path = sys.argv[1]
    out_path = sys.argv[2]

    data_occ = {}
    for record in SeqIO.parse(fq_path, "fastq"):
        idx, ab = record.id.strip('\n').split('#')
        ab = int(ab)
        l = len(record)
        data_occ[(l,ab)] = data_occ[(l,ab)]+1 if (l,ab) in data_occ else 1

    data = {"Length":[], "Abundance":[]}
    for (l,ab),occ in data_occ.items():
        # +1 because log(1) = 0
        data["Length"] += [l]*(occ+1)
        data["Abundance"] += [ab]*(occ+1)

    df = pd.DataFrame(data = data)
    p = sns.jointplot(data=data, x="Length", y="Abundance", color="seagreen",
                      marginal_ticks = True,
                      marginal_kws = dict(bins=8000, fill=False))
    p.ax_marg_x.yaxis.label.set_visible(True)
    p.ax_marg_y.xaxis.label.set_visible(True)
    p.ax_marg_x.set_yscale("log")
    p.ax_marg_x.set_ylabel("Count (log)")
    p.ax_marg_y.set_xscale("log")
    p.ax_marg_y.set_xlabel("Count (log)")

    axins = p.ax_joint.inset_axes([0.6, 0.44, 0.39, 0.55])
    minL = []
    minA = []
    for i in range(0,len(data["Length"])):
        if data["Length"][i] <= 500:
            minL.append(data["Length"][i])
            minA.append(data["Abundance"][i])
    axins.scatter(minL, minA, color="seagreen", alpha=1, linewidths = 1, edgecolors="white")
    p.ax_joint.indicate_inset_zoom(axins, edgecolor="grey", ls="--")
    plt.savefig(out_path)
    # plt.show()

def plot_nal_by_err():
    bampath1 = sys.argv[1]
    bampath2 = sys.argv[2]
    fq_path = sys.argv[3]
    out_path = sys.argv[4]

    tot_strings = 0
    for record in SeqIO.parse(fq_path, "fastq"):
        tot_strings += 1
    print("Tot strings: ", tot_strings)

    data1 = []
    bamfile1 = pysam.AlignmentFile(bampath1, "rb")
    for al in bamfile1.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        data1.append(bad_bases)
    print("bam1 limits: ", min(data1), max(data1))

    data2 = []
    bamfile2 = pysam.AlignmentFile(bampath2, "rb")
    for al in bamfile2.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        data2.append(bad_bases)
    print("bam2 limits: ", min(data2), max(data2))

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
    plt.savefig(out_path)
    # plt.show()

def plot_covering():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]
    out_path = sys.argv[3]

    covering = {}
    not_covering = {}
    
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, iden, count = line[3], int(line[4]), int(line[-1])
        if count == 0:
            not_covering[idx] = iden
        else:
            covering[idx] = iden

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, iden, count = line[3], int(line[4]), int(line[-1])
        if count == 0:
            if idx in covering:
                pass
            elif idx in not_covering:
                not_covering[idx] = max(iden, not_covering[idx])
            else:
                not_covering[idx] = iden
        else:
            if idx in covering:
                covering[idx] = max(iden, covering[idx])
            elif idx in not_covering:
                not_covering.pop(idx)
                covering[idx] = iden
            else:
                covering[idx] = iden

    fig, ax1 = plt.subplots(1, 1, tight_layout=True)
    nbins = 7
    ax1.hist([list(covering.values()), list(not_covering.values())],
             bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
             histtype="bar", cumulative=False,
             color=["seagreen", "darkorchid"], label=["Covering", "Not covering"])

    ax1.set_xlabel("# Errors")
    ax1.set_ylabel("# Primary Alignments")

    plt.legend(loc=1)
    plt.savefig(out_path)
    # plt.show()

if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "lsamples":
        plot_samplelengths()
    elif mode == "corr":
        plot_correlation()
    elif mode == "albyerr":
        plot_nal_by_err()
    elif mode == "covplot":
        plot_covering()
