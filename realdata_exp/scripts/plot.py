import sys

from Bio import SeqIO
import pysam

import numpy as np
import pandas as pd
import matplotlib
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
        if data["Abundance"][i] >= 5 and data["Length"][i] <= 500:
            minL.append(data["Length"][i])
            minA.append(data["Abundance"][i])
    axins.scatter(minL, minA, color="seagreen", alpha=1, linewidths = 1, edgecolors="white")
    axins.set_xlim(0, 500)
    axins.set_ylim(0, max(data["Abundance"])+100)
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
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        # bbmap:
        good_bases = al.get_cigar_stats()[0][7]
        # minimap2:
        # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data1.append(errors)
    print("bam1 limits: ", min(data1), max(data1))

    data2 = []
    bamfile2 = pysam.AlignmentFile(bampath2, "rb")
    for al in bamfile2.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        # bbmap:
        good_bases = al.get_cigar_stats()[0][7]
        # minimap2:
        # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data2.append(errors)
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
    plt.xlabel("# Base Differences")

    ax2.legend(loc=4, fancybox=True, shadow=True)
    plt.savefig(out_path)
    # plt.show()

def plot_covering():
    bed1_path = sys.argv[1]
    bed2_path = sys.argv[2]
    sam1_path = sys.argv[3]
    sam2_path = sys.argv[4]
    out_path = sys.argv[5]

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
        # bbmap:
        good_bases = al.get_cigar_stats()[0][7]
        # minimap2:
        # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data1.append(errors)

    data2 = []
    bamfile2 = pysam.AlignmentFile(sam2_path, "rb")
    for al in bamfile2.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        # bbmap:
        good_bases = al.get_cigar_stats()[0][7]
        # minimap2:
        # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
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

    plt.savefig(out_path)

# def plot_covering():
#     bed_path_1 = sys.argv[1]
#     bed_path_2 = sys.argv[2]
#     out_path = sys.argv[3]

#     alignments = {}
#     for line in open(bed_path_1):
#         line = line.strip('\n').split('\t')
#         idx, err, count = line[3], int(line[4]), int(line[-1])
#         if idx not in alignments:
#             alignments[idx] = (False, float("inf"))
#         if count > 0:
#             alignments[idx] = (True, min(alignments[idx][1], err))
#         else:
#             if not alignments[idx][0]:
#                 alignments[idx] = (False, min(alignments[idx][1], err))

#     for line in open(bed_path_2):
#         line = line.strip('\n').split('\t')
#         idx, err, count = line[3], int(line[4]), int(line[-1])
#         if idx not in alignments:
#             alignments[idx] = (False, float("inf"))
#         if count > 0:
#             alignments[idx] = (True, min(alignments[idx][1], err))
#         else:
#             if not alignments[idx][0]:
#                 alignments[idx] = (False, min(alignments[idx][1], err))

#     covering = []
#     not_covering = []
#     for idx, (cov_flag, err) in alignments.items():
#         if cov_flag:
#             covering.append(err)
#         else:
#             not_covering.append(err)

#     print("Covering:", len(covering))
#     print("Limits:", min(covering), max(covering))
#     print("---")
#     print("Not covering:", len(not_covering))
#     print("Limits:", min(not_covering), max(not_covering))

#     fig, ax1 = plt.subplots(1, 1, tight_layout=True)
#     nbins = 7
#     ax1.hist([covering, not_covering],
#              bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
#              histtype="bar", cumulative=False,
#              color=["seagreen", "darkorchid"], label=["Covering", "Not covering"])

#     ax1.set_xlabel("# Base Differences")
#     ax1.set_ylabel("# Primary Alignments")

#     plt.legend(loc=1)
#     plt.savefig(out_path)
#     # plt.show()

def plot_covering_abundances():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]
    out_path = sys.argv[3]

    covering = set()
    not_covering = set()
    
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if count == 0:
            if idx in covering:
                pass
            else:
                not_covering.add(idx)
        else:
            if idx in not_covering:
                not_covering.remove(idx)
            covering.add(idx)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if count == 0:
            if idx in covering:
                pass
            elif idx in not_covering:
                pass
            else:
                not_covering.add(idx)
        else:
            if idx in covering:
                pass
            elif idx in not_covering:
                not_covering.remove(idx)
                covering.add(idx)
            else:
                covering.add(idx)

    print("Covering:", len(covering))
    print("Non-Covering:", len(not_covering))
    print("Total:", len(covering) + len(not_covering))

    data = {"Type" : ["Covering" for _ in covering] + ["Non-Covering" for _ in not_covering],
            "Abundance" : [int(idx.split('#')[1]) for idx in covering] + [int(idx.split('#')[1]) for idx in not_covering]}
    df = pd.DataFrame(data)

    sns.histplot(data=df, x="Abundance", hue="Type", element="step",
                 binrange=(0,50),
                 bins=range(0,51))

    # plt.savefig(out_path)
    plt.show()

# def plot_pmatches_len():
#     bampath = sys.argv[1]
#     out_path = sys.argv[2]

#     data = []
#     bamfile = pysam.AlignmentFile(bampath, "rb")
#     for al in bamfile.fetch():
#         good_bases = al.get_cigar_stats()[0][7]
#         bad_bases = al.query_length - good_bases
#         deletions = al.get_cigar_stats()[0][2]
#         errors = bad_bases + deletions
#         if errors == 0:
#             data.append(al.query_length)
#     print("Limits: ", min(data), max(data))

#     nbins = 501
#     fig, ax1 = plt.subplots(1, 1, sharey = True)

#     ax1.hist(data, bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left", color="seagreen")
#     ax1.set_xlabel("Length")
#     ax1.set_ylabel("# specific strings")
#     plt.savefig(out_path)

def recbylen():
    fpath = sys.argv[1]
    out_path = sys.argv[2]

    data = {}
    for line in open(fpath):
        idx, l, ov = line.strip('\n').split(' ')
        l, ov = int(l), int(ov)
        if l not in data:
            data[l] = [0,0]
        data[l][0] += ov > 0
        data[l][1] += 1

    Xs = sorted(data.keys())
    Ys1 = [data[x][0]/data[x][1] for x in Xs]

    total = []
    for l,(_,tot) in data.items():
        total += [l]*(tot+1)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, tight_layout=True) #, figsize=(13, 7))

    sns.scatterplot(x = Xs, y = Ys1, hue = Ys1, ax = ax1, legend = False,
                    palette = sns.dark_palette("seagreen", as_cmap=True), linewidth=0, alpha = 1, s=10)

    sns.histplot(total, bins=250, log_scale = (False,True), color = "seagreen", element="poly")
    ax2.invert_yaxis()

    ax1.set_ylabel("Recall")
    ax2.set_xlabel("Variation size")
    ax2.set_ylabel("Count (log)")

    plt.savefig(out_path)
    # plt.show()

def recbylen2():
    fpath = sys.argv[1]
    out_path = sys.argv[2]

    data = {}
    for line in open(fpath):
        idx, l, ov = line.strip('\n').split(' ')
        l, ov = int(l), int(ov)
        if l not in data:
            data[l] = [0,0]
        data[l][0] += ov > 0
        data[l][1] += 1

    Xs = sorted(data.keys())
    Ys1 = [data[x][0]/data[x][1] for x in Xs]

    found = []
    for l,(f,_) in data.items():
        found += [l]*(f+1)
    
    total = []
    for l,(_,tot) in data.items():
        total += [l]*(tot+1)

    p1 = sns.histplot(total, bins=500, log_scale = (False,True), color = "black", element="poly", label="Not found")
    p2 = sns.histplot(found, bins=500, log_scale = (False,True), color = "lime", alpha=0.75, element="poly", label="Found")

    p2.set_xlabel("Variation size")
    p2.set_ylabel("Count (log)")

    plt.legend()
    # plt.savefig(out_path)
    plt.show()

def repmask():
    font = {"size" : 13}
    matplotlib.rc("font", **font)
    
    fa_path = sys.argv[1]
    rm_path = sys.argv[2]
    out_path = sys.argv[3]
    min_l = int(sys.argv[4]) if len(sys.argv) == 5 else 0

    lens = {}
    L = 0
    N = 0
    for record in SeqIO.parse(fa_path, "fasta"):
        lens[record.id] = len(record)
        seq = str(record.seq).upper()
        l = len(seq)
        if l >= min_l:
            n = seq.count('N')
            L += l
            N += n
    print(f"Masked: {N}/{L} ({N/L})")

    labels = ["SINE", "LINE", "Satellite", "Simple Repeat", "Low Complex.", "Other", "None"]
    rtypes = {"SINE" : "SINE",
              "LINE" : "LINE",
              "Satellite" : "Satellite",
              "Simple_repeat" : "Simple Repeat",
              "Low_complexity" : "Low Complex.",
              "Other" : "Other"}
    data = {"Repeat Type" : [], "Specific String Length" : []}
    used_ridx = set()
    for line in open(rm_path, 'r'):
        line = line.split()
        # print("\t".join(line))
        if len(line) == 0 or line[0] == "SW" or line[0] == "score":
            continue
        ridx, rtype = line[4], line[10].split('/')[0]
        if rtype not in rtypes:
            rtype = "Other"
        if lens[ridx] >= min_l:
            data["Repeat Type"].append(rtypes[rtype])
            data["Specific String Length"].append(lens[ridx])
            used_ridx.add(ridx)

    for ridx in lens:
        if ridx not in used_ridx:
            if lens[ridx] >= min_l:
                data["Repeat Type"].append("None")
                data["Specific String Length"].append(lens[ridx])

    df = pd.DataFrame(data)

    for t in list(rtypes.values()) + ["None"]:
        print(t, len(df.loc[df["Repeat Type"] == t])/len(df)*100)

    dplot = sns.histplot(data=df, x="Specific String Length", hue="Repeat Type", multiple="stack", bins=30, hue_order=labels)

    plt.tight_layout()
    plt.savefig(out_path)
    # plt.show()

def plot_nalignment():
    sampath = sys.argv[1]
    out_path = sys.argv[2]

    aligns = {}
    sam = pysam.AlignmentFile(sampath, "r")
    for al in sam.fetch():
        aligns[al.query_name] = aligns[al.query_name]+1 if al.query_name in aligns else 1

    data = {}
    for v in aligns.values():
        if v not in data:
            data[v] = 0
        data[v] += 1

    Xs = sorted(data.keys())
    Ys = [data[x] for x in Xs]

    fig, ax = plt.subplots()
    ax.bar(Xs, Ys)
    ax.set_title("Alignments to HG00733 contigs")
    ax.set_xlabel("# alignments")
    ax.set_ylabel("Count")
    plt.savefig(out_path)

if __name__ == "__main__":
    mode = sys.argv.pop(1) if len(sys.argv) > 2 else ""
    if mode == "lsamples":
        plot_samplelengths()
    elif mode == "corr":
        plot_correlation()
    elif mode == "albyerr":
        plot_nal_by_err()
    elif mode == "cov":
        plot_covering()
    # elif mode == "pmatches":
    #     plot_pmatches_len()
    elif mode == "recbylen":
        recbylen()
    elif mode == "repmask":
        repmask()
    elif mode == "nal":
        plot_nalignment()
    elif mode == "covab":
        plot_covering_abundances()
    else:
        print("python3 plot.py")
        print("\t\tlsamples <fq1> <fq2> <out>")
        print("\t\tcorr <fq> <out>")
        print("\t\talbyerr <bam1> <bam2> <fq> <out>")
        print("\t\tcov <bed1> <bed2> <out>")
        print("\t\tpmatches <bam> <out>")
        print("\t\trecbylen <recall.list> <out>")
        print("\t\trepmask <fa> <rm.out> <out> [min_l]")
        print("\t\tnal <bam> <out>")
