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

def main():
    f = open(sys.argv[1]) if len(sys.argv)>1 else sys.stdin
    eprint(f"Parsing...")
    X = []
    for line in f:
        x = int(line.strip('\n'))
        if x>1:
            X.append(x)
    f.close()

    eprint("Computing statistics...")
    X = np.array(X)
    df = pd.Series(X)
    print(df.describe())

    eprint("Plotting...")
    bins_list = list(range(max(X)+1))
    fig, (ax1) = plt.subplots(1,1)
    sns.distplot(X,
                 color = "seagreen",
                 bins=bins_list,
                 kde = False,
                 norm_hist = False,
                 ax = ax1)
    eprint("Saving dist plot...")
    plt.savefig("simple_hist.pdf")
    eprint("Done.")

if __name__ == "__main__":
    main()
