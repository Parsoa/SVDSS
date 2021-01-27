import sys

from Bio import SeqIO
import pysam

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_sv_lengths():
    bed_path1 = sys.argv[1]
    out_path = sys.argv[2]

    data1 = []
    for line in open(bed_path1).readlines():
        if line[0] != '#':
            line = line.split()
            length = int(line[4])
            data1.append(length)

    data = {"SVs" : ["INS"] * len(data1) + ["DEL"] * len(data2),
            "Lengths" : data1 + data2 }
    df = pd.DataFrame(data = data)

    plt.figure(figsize=(24, 8))
    dplot = sns.displot(df, x = "Lengths", hue = "SVs", fill = True, legend = False)
    plt.savefig(out_path)

def plot_abundance():
    bed_path1 = sys.argv[1]
    bed_path2 = sys.argv[2]
    out_path = sys.argv[3]

    data1 = []
    for line in open(bed_path1).readlines():
        line = line.split()
        abundance = int(line[4])
        if abundance < 50:
            data1.append(abundance)
    data2 = []
    for line in open(bed_path2).readlines():
        line = line.split()
        abundance = int(line[4])
        if abundance < 50:
            data2.append(abundance)

    #data2 = [s / len(data2) for s in data2]
    #data1 = [s / len(data1) for s in data1]

    #data = {"Strings" : ["SV"]*len(data1) + ["non-SV"]*len(data2),
    #"Abundances" : data1 + data2}
    data = {"Strings" : ["Long"] * len(data2),
            "Abundance" : data2}
    df = pd.DataFrame(data = data)

    nbins = 15
    dplot = sns.histplot(df, x="Abundance", hue="Strings", fill=True, legend=False, bins=np.arange(0, nbins+1))
    #dplot.ax.legend(labels=["non-SV", "SV"], loc=2)

    plt.savefig('plots/' + out_path)
    # plt.show()

def plot_unmapped_length_both():
    bed_path1 = sys.argv[1]
    bed_path2 = sys.argv[2]
    out_path = sys.argv[3]

    data1 = []
    for line in open(bed_path1).readlines():
        line = line.split()
        data1.append(len(line[3]))
    data2 = []
    for line in open(bed_path2).readlines():
        line = line.split()
        data2.append(len(line[3]))

    #data1 = [max(d, 500) for d in data1]
    #data2 = [max(d, 500) for d in data2]

    data = {"Strings" : ["Short"]*len(data1) + ["Long"]*len(data2),
            "Lengths" : data1 + data2}
    df = pd.DataFrame(data = data)

    dplot = sns.histplot(df, x="Lengths", hue="Strings", fill=True, legend=False, binwidth = 10)
    plt.ylim(0, 500)

    plt.savefig(out_path)
    # plt.show()

def plot_unmapped_length():
    bed_path1 = sys.argv[1]
    out_path = sys.argv[2]

    data1 = []
    for line in open(bed_path1).readlines():
        line = line.split()
        data1.append(len(line[3]))

    df = pd.DataFrame(data = data1)
    dplot = sns.displot(df, legend = False)
    plt.savefig(out_path)

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
        if data["Abundance"][i] >= 5 and data["Length"][i] <= 500:
            minL.append(data["Length"][i])
            minA.append(data["Abundance"][i])
    axins.scatter(minL, minA, color="seagreen", alpha=1, linewidths = 1, edgecolors="white")
    axins.set_xlim(0, 500)
    axins.set_ylim(0, max(data["Abundance"])+100)
    p.ax_joint.indicate_inset_zoom(axins, edgecolor="grey", ls="--")
    plt.savefig(out_path)
    # plt.show()

def plot_readlength():
    fq_path = sys.argv[1]
    out_path = sys.argv[2]
    data = []
    for record in SeqIO.parse(fq_path, "fastq"):
        data.append(len(record))
        if len(data) > 100000:
            break
    df = pd.DataFrame(data = data)
    dplot = sns.histplot(data = df)
    plt.savefig(out_path)

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

    ax1.set_xlabel("# Base Differences")
    ax1.set_ylabel("# Primary Alignments")

    plt.legend(loc=1)
    plt.savefig(out_path)

def plot_nal_by_err():
    fastapth = sys.argv[1]
    out_path = sys.argv[2]

    tot_strings = len(open(fastapth).readlines()) / 2
    print("Tot strings: ", tot_strings)

    data1 = []
    for bampath1 in ["short.sorted.bam"]:
    #for bampath1 in [u"short.sorted.bam", "long.bam"]:
        bamfile1 = pysam.AlignmentFile(bampath1, "rb")
        for al in bamfile1.fetch():
            if al.is_secondary or al.is_unmapped or al.is_supplementary:
                continue
            good_bases = al.get_cigar_stats()[0][7]
            bad_bases = al.query_length - good_bases
            deletions = al.get_cigar_stats()[0][2]
            errors = bad_bases + deletions
            data1.append(errors)
    print("bam1 limits: ", min(data1), max(data1))
    print(len(data1))

    data2 = []
    for bampath2 in ["short.mother.sorted.bam"]:
    #for bampath2 in ["long.mother.bam"]:
        bamfile2 = pysam.AlignmentFile(bampath2, "rb")
        for al in bamfile2.fetch():
            if al.is_secondary or al.is_unmapped or al.is_supplementary:
                continue
            good_bases = al.get_cigar_stats()[0][7]
            bad_bases = al.query_length - good_bases
            deletions = al.get_cigar_stats()[0][2]
            errors = bad_bases + deletions
            data2.append(errors)
    print("bam2 limits: ", min(data2), max(data2))
    print(len(data2))

    data3 = []
    for bampath3 in ["short.father.sorted.bam"]:
    #for bampath3 in ["long.father.bam"]:
        bamfile3 = pysam.AlignmentFile(bampath3, "rb")
        for al in bamfile3.fetch():
            if al.is_secondary or al.is_unmapped or al.is_supplementary:
                continue
            good_bases = al.get_cigar_stats()[0][7]
            bad_bases = al.query_length - good_bases
            deletions = al.get_cigar_stats()[0][2]
            errors = bad_bases + deletions
            data3.append(errors)
    print("bam3 limits: ", min(data3), max(data3))
    print(len(data3))

    nbins = 6
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey = True, tight_layout=True, figsize=(13, 5))
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey = True, tight_layout=True, figsize=(8, 4))

    ax1.axhline(tot_strings, 0, 100, linewidth=1, color='r', label="Total specific")
    ax1.hist(data1, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", alpha=0.75, cumulative = True, histtype="step", label="Cumulative count")
    ax1.hist(data1, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", histtype="stepfilled", label="Count")
    ax1.set_title("(c) Child", fontsize = 13)
    ax1.set_ylabel("# Primary Alignments", fontsize = 13)

    ax2.axhline(tot_strings, 0, 100, linewidth=1, color='r', label="Total specific")
    ax2.hist(data2, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", alpha=0.75, cumulative = True, histtype="step", label="Cumulative count")
    ax2.hist(data2, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", histtype="stepfilled", label="Count")
    ax2.set_title("(d) Mother", fontsize = 13)

    #ax3.axhline(tot_strings, 0, 100, linewidth=1, color='r', label="Total specific")
    #ax3.hist(data3, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", alpha=0.75, cumulative = True, histtype="step", label="Cumulative count")
    #ax3.hist(data3, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen", histtype="stepfilled", label="Count")
    #ax3.set_title("Father")
    
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xticks()
    plt.xlabel("# Base Differences", fontsize = 13)

    ax1.legend(loc=4, fancybox=True, shadow=True)
    plt.savefig(out_path, dpi = 300)

def repmask():
    fa_path = sys.argv[1]
    rm_path = sys.argv[2]
    out_path = sys.argv[3]
    min_l = int(sys.argv[4]) if len(sys.argv) == 4 else 0

    lens = {}
    for record in SeqIO.parse(fa_path, "fasta"):
        lens[record.id] = len(record)

    rtypes = ["SINE", "LINE", "Satellite", "Simple_repeat", "Low_complexity", "Other"]
    data = {"rtype":[], "len":[]}
    used_ridx = set()
    for line in open(rm_path, 'r'):
        line = line.split()
        if len(line) == 0 or line[0] == "SW" or line[0] == "score":
            continue
        ridx, rtype = line[4], line[10].split('/')[0]
        if rtype not in rtypes:
            rtype = "Other"
        if ridx in lens:
            if lens[ridx] >= min_l:
                data["rtype"].append(rtype)
                data["len"].append(lens[ridx])
                used_ridx.add(ridx)

    for ridx in lens:
        if ridx not in used_ridx:
            if lens[ridx] >= min_l:
                data["rtype"].append("None")
                data["len"].append(lens[ridx])

    df = pd.DataFrame(data)
    dplot = sns.histplot(data=df, x="len", hue="rtype", multiple="stack")
    plt.savefig('plots/' + out_path)

def plot_covering():

    # plot (hardcoded) parameters
    nbins = 7
    barw = 0.35
    Xs = list(range(0,nbins))
    
    # Parsing BEDs
    alignments = {}
    for line in open(bed1_path):
        line = line.strip('\n').split('\t')
        idx, err, count = line[3], int(line[4]), int(line[-1])
        if idx not in alignments:
            alignments[idx] = (False, float("inf"))
        if count > 0:
            alignments[idx] = (True, min(alignments[idx][1], err))
        else:
            if not alignments[idx][0]:
                alignments[idx] = (False, min(alignments[idx][1], err))

    for line in open(bed2_path):
        line = line.strip('\n').split('\t')
        idx, err, count = line[3], int(line[4]), int(line[-1])
        if idx not in alignments:
            alignments[idx] = (False, float("inf"))
        if count > 0:
            alignments[idx] = (True, min(alignments[idx][1], err))
        else:
            if not alignments[idx][0]:
                alignments[idx] = (False, min(alignments[idx][1], err))

    covering = []
    not_covering = []
    for idx, (cov_flag, err) in alignments.items():
        if cov_flag:
            covering.append(err)
        else:
            not_covering.append(err)

    print("Covering:", len(covering))
    print("Limits:", min(covering), max(covering))
    print("---")
    print("Not covering:", len(not_covering))
    print("Limits:", min(not_covering), max(not_covering))

    # Parsing BAMs
    data1 = []
    bamfile1 = pysam.AlignmentFile(sam1_path, "rb")
    for al in bamfile1.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data1.append(errors)

    data2 = []
    bamfile2 = pysam.AlignmentFile(sam2_path, "rb")
    for al in bamfile2.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data2.append(errors)

    data1_dict = {}
    for e in data1:
        if e < nbins:
            data1_dict[e] = data1_dict[e]+1 if e in data1_dict else 1
    data2_dict = {}
    for e in data2:
        if e < nbins:
            data2_dict[e] = data2_dict[e]+1 if e in data2_dict else 1

    print("bam1 limits: ", min(data1), max(data1))
    print("bam2 limits: ", min(data2), max(data2))

    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey = True, tight_layout=True, figsize=(13, 7))

    ax1.hist([covering, not_covering],
             bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
             histtype="bar", cumulative=False,
             color=["steelblue", "darkorchid"], label=["Covering", "Not covering"])
    ax1.set_ylabel("# Primary Alignments")
    ax1.legend()

    ax2.bar(Xs, [data1_dict[x] for x in Xs], color="plum", width=-barw, align="edge", label="HG00733")
    ax2.bar(Xs, [data2_dict[x] for x in Xs], color="pink", width=barw, align="edge", label="NA19240")
    ax2.legend()

    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
    plt.xticks()
    plt.xlabel("# Base Differences")

    plt.savefig(out_path, dpi = 300)

def plot_cutoff():
    events = 17595
    data_5x = {'cutoff': [2, 3, 4, 5, 6],
            'mapped': [1698507, 1681697, 1670331, 1649708, 1606824],
            'unmapped': [2181981, 539619, 135258, 41708, 14624], 
            'covered': [17354, 17344, 17291, 17190, 16944]}
    data_10x = {'cutoff': [2, 3, 4, 5, 6],
            'mapped': [1698507, 1681697, 1670331, 1649708, 1606824],
            'unmapped': [2181981, 539619, 135258, 41708, 14624], 
            'covered': [17354, 17344, 17291, 17190, 16944]}
    data_20x = {'cutoff': [2, 3, 4, 5, 6],
            'mapped': [1745333, 1712197, 1702999, 1700609, 1699879],
            'unmapped': [6245298, 1957703, 486772, 133528, 43215], 
            'covered': [17414, 17414, 17414, 17413, 17413]}
    data_30x = {'cutoff': [2, 3, 4, 5, 6],
            'mapped': [1767781, 1714017, 1695793, 1690675, 1689266],
            'unmapped': [12136901, 4607079, 1286537, 361469, 131090], 
            'covered': [17371, 17369, 17369, 17368, 17368]}
    # Nono-overlapping
    #data_5x  = {'cutoff': [2, 3, 4, 5, 6], 'precision': [0.81, 0.97, 0.99, 0.99, 0.99], 'recall': [0.94, 0.82, 0.65, 0.48, 0.34]}
    #data_10x = {'cutoff': [2, 3, 4, 5, 6], 'precision': [0.66, 0.94, 0.99, 0.99, 0.99], 'recall': [0.98, 0.97, 0.95, 0.89, 0.78]}
    #data_20x = {'cutoff': [2, 3, 4, 5, 6], 'precision': [0.39, 0.82, 0.97, 0.99, 0.99], 'recall': [0.98, 0.98, 0.98, 0.98, 0.98]}
    #data_30x = {'cutoff': [2, 3, 4, 5, 6], 'precision': [0.24, 0.64, 0.91, 0.97, 0.97], 'recall': [0.99, 0.99, 0.99, 0.99, 0.99]}

    data = {'cutoff': data_5x['cutoff'] + data_10x['cutoff'] + data_20x['cutoff'] + data_30x['cutoff'],
            'precision': [d['mapped'][i] / (d['mapped'][i] + d['unmapped'][i]) for d in [data_5x, data_10x, data_20x, data_30x] for i in [0, 1, 2, 3, 4]],
            'recall': [d['covered'][i] / events for d in [data_5x, data_10x, data_20x, data_30x] for i in [0, 1, 2, 3, 4]],
            'coverage': [5] * 5 + [10] * 5 + [20] * 5 + [30] * 5}
    
    #data = {'cutoff': data_5x['cutoff'] + data_10x['cutoff'] + data_20x['cutoff'] + data_30x['cutoff'],
    #        'precision': data_5x['precision'] + data_10x['precision'] + data_20x['precision'] + data_30x['precision'],
    #        'recall': data_5x['recall'] + data_10x['recall'] + data_20x['recall'] + data_30x['recall'],
    #        'coverage': [5] * 5 + [10] * 5 + [20] * 5 + [30] * 5}

    df = pd.DataFrame(data=data)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey = True, tight_layout=True, figsize=(8, 4))
    sns.lineplot(ax = ax1, data=data, x="cutoff", y="precision", hue="coverage", palette="Dark2", legend = False)
    ax1.set_title("(a) Precision", fontsize = 13)
    ax1.set_ylabel("")
    ax1.set_xlabel("")
    ax1.set_xticks(list(range(2, 7)))
    #plt.savefig('precision.png')
    sns.lineplot(ax = ax2, data=data, x="cutoff", y="recall", hue="coverage", palette="Dark2")
    ax2.set_title("(b) Recall", fontsize = 13)
    ax2.set_ylabel("")
    ax2.set_xlabel("")
    ax2.set_xticks(list(range(2, 7)))
    
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
    plt.xticks()
    plt.xlabel('Cutoff', fontsize = 13)
    plt.savefig('cutoff.png', dpi = 300)

if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "svlen":
        plot_sv_lengths()
    if mode == "abund":
        plot_abundance()
    elif mode == "corr":
        plot_correlation()
    elif mode == "covplot":
        plot_covering()
    elif mode == "repmask":
        repmask()
    elif mode == "nalbyerr":
        plot_nal_by_err()
    elif mode == "readlength":
        plot_readlength()
    elif mode == "length_unmapped":
        plot_unmapped_length()
    elif mode == "length_unmapped_both":
        plot_unmapped_length_both()
    elif mode == "cutoff":
        plot_cutoff()
