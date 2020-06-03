import sys, os, glob
import time

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

import seaborn as sns

def eprint(s):
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S %m/%d", t)
    print(f"[{current_time}] {s}", file=sys.stderr)

def parse_file(fpath):
    eprint(f"Parsing file {fpath}...")
    X = []
    for line in open(fpath):
        x = int(line.strip('\n'))
        X.append(x)
    return np.array(X)

def main():
    fpath1 = sys.argv[1]
    fpath2 = sys.argv[2]
    out_prefix = sys.argv[3] if len(sys.argv) >= 4 else "_out"

    len1 = parse_file(fpath1)
    len2 = parse_file(fpath2)
    len2_masked = [x for x in len2 if x>=31]
    len2_masked_2 = [x for x in len2 if x>=2500]

    eprint("Computing statistics...")
    df = pd.Series(len2)
    print(df.describe())

    eprint(f"Plotting...")
    bins_list = list(range(0, 25000))
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, sharex=True)
    ax1.hist(len1,
             color = "darkorchid",
             bins=bins_list)
    ax2.hist(len2_masked,
             color = "seagreen",
             bins=bins_list)
    ax3.hist(len2_masked_2,
             color = "seagreen",
             bins=bins_list)
    ax4.hist(len2,
             color = "steelblue",
             log=True,
             bins=bins_list)
    ax4.set_xlabel("Read length")

    eprint("Saving hist plot...")
    plt.savefig(f"{out_prefix}.hist.pdf")

    eprint("Done.")

def main_multi():
    fq_path = sys.argv[1]
    fa_paths = sys.argv[2]

    eprint(f"Parsing sample {fq_path}...")
    read_len = []
    for line in open(fq_path):
        read_len.append(int(line.strip('\n')))

    cut_len = {}
    for fa_path in glob.glob(fa_paths):
        eprint(f"Parsing sample {fa_path}...")
        m = int(fa_path.split('.')[-4][1:])
        cut_len[m] = []
        for line in open(fa_path):
            cut_len[m].append(int(line.strip('\n')))


    read_len = np.array(read_len)
    bins_list = list(range(0, 25000))
    # fig, axes = plt.subplots(2,2)
    # axes = axes.reshape(-1)
    colors = ["lime", "aqua", "violet", "coral", "lightpink"]
    i = 0
    for k in [5]: #[5,10,25,50]:
        data = np.array(cut_len[k])

        eprint(f"Plotting {k}...")
        fig, (ax1,ax2) = plt.subplots(2,1)
        sns.distplot(data,
                     color = "black",
                     bins=bins_list,
                     ax = ax1)
        sns.distplot(data,
                     color = colors[i],
                     bins=bins_list,
                     ax = ax2)
        i+=1
        eprint("Saving plot...")
        plt.savefig(f"hist.{k}.pdf")
    eprint("Done.")

if __name__ == "__main__":
    main()
