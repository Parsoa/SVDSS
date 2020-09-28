import sys

from Bio import SeqIO

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_lengths():
    # each line is the length of a read
    # eg sed -n '2~4p' <fastq> | awk '{print length}'
    fpath1 = sys.argv[1]
    fpath2 = sys.argv[2]

    data1 = []
    for line in open(fpath1):
        data1.append(int(line))

    data2 = []
    for line in open(fpath2):
        data2.append(int(line))

    data = {"Sample" : ["HG00733"]*len(data1) + ["NA19240"]*len(data2),
            "Lengths" : data1 + data2}
    df = pd.DataFrame(data = data)
    
    sns.displot(df, x="Lengths", hue="Sample", fill=True)
    plt.savefig("samples.len.png")
    # plt.show()


def plot_correlation():
    fq_path = sys.argv[1]

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
    plt.savefig("corr.zoom.png")
    # plt.show()

if __name__ == "__main__":
    plot_correlation()
