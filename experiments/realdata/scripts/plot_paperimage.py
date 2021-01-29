import sys

from Bio import SeqIO
import pysam

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_al_errors(al, isbbmap):
    if isbbmap:
        good_bases = al.get_cigar_stats()[0][7]
    else: # minimap2:
        good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
    bad_bases = al.query_length - good_bases
    deletions = al.get_cigar_stats()[0][2]
    return bad_bases + deletions

def extract_from_bam(bampath, isbbmap, cutoff = float("inf")):
    data = []
    bamfile = pysam.AlignmentFile(bampath, "rb")
    for al in bamfile.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        errors = get_al_errors(al, isbbmap)
        if errors < cutoff:
            data.append(errors)
    return data

def main():
    nbins = 6

    fq_path = sys.argv[1]
    hgcontigsshort_sam = sys.argv[2]
    hgcontigslong_sam = sys.argv[3]
    nacontigsshort_sam = sys.argv[4]
    nacontigslong_sam = sys.argv[5]
    shortuncovering_txt = sys.argv[6]
    longuncovering_txt = sys.argv[7]
    hghaploshort_sam = sys.argv[8]
    hghaplolong_sam = sys.argv[9]
    subhgcontigsshort_sam = sys.argv[10]
    subhgcontigslong_sam = sys.argv[11]
    subnacontigsshort_sam = sys.argv[12]
    subnacontigslong_sam = sys.argv[13]
    out_path = sys.argv[14]

    # Parsing sample
    tot_strings = 0
    for line in open(fq_path):
        tot_strings += 1
    tot_strings /= 4 # Assuming FASTQ
    print("Tot strings: ", tot_strings)

    # Parsing alignments to contigs
    hgcontigsshort_nalbyerr = extract_from_bam(hgcontigsshort_sam, True, nbins)
    hgcontigslong_nalbyerr = extract_from_bam(hgcontigslong_sam, False, nbins)
    hgcontigs_nalbyerr = hgcontigsshort_nalbyerr + hgcontigslong_nalbyerr
    print("On HG contigs:", len(hgcontigs_nalbyerr))

    nacontigsshort_nalbyerr = extract_from_bam(nacontigsshort_sam, True, nbins)
    nacontigslong_nalbyerr = extract_from_bam(nacontigslong_sam, False, nbins)
    nacontigs_nalbyerr = nacontigsshort_nalbyerr + nacontigslong_nalbyerr
    print("On NA contigs:", len(nacontigs_nalbyerr))

    # Parsing non haplo-compatible strings (after alignment to HG haplo)
    nonhaplocomp_idxs = set()
    for line in open(shortuncovering_txt):
        nonhaplocomp_idxs.add(line.strip('\n'))
    for line in open(longuncovering_txt):
        nonhaplocomp_idxs.add(line.strip('\n'))
    print("Non Haplo-Compatible:", len(nonhaplocomp_idxs))

    # Parsing haplo sam
    haplocomp = []
    nonhaplocomp = []
    for al in pysam.AlignmentFile(hghaploshort_sam, 'r').fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        errors = get_al_errors(al, True)
        if errors < nbins:
            if al.query_name in nonhaplocomp_idxs:
                nonhaplocomp.append(errors)
            else:
                haplocomp.append(errors)
    for al in pysam.AlignmentFile(hghaplolong_sam, 'r').fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        errors = get_al_errors(al, False)
        if errors < nbins:
            if al.query_name in nonhaplocomp_idxs:
                nonhaplocomp.append(errors)
            else:
                haplocomp.append(errors)
    print("Haplo-Aware:", len(haplocomp))
    print("Non Haplo-Aware:", len(nonhaplocomp))

    # Parsing non haplo-compatible on contigs
    subhgcontigsshort_nalbyerr = extract_from_bam(subhgcontigsshort_sam, True, nbins)
    subhgcontigslong_nalbyerr = extract_from_bam(subhgcontigslong_sam, False, nbins)
    subhgcontigs_nalbyerr = subhgcontigsshort_nalbyerr + subhgcontigslong_nalbyerr
    print("On HG contigs (sub):", len(subhgcontigs_nalbyerr))

    subnacontigsshort_nalbyerr = extract_from_bam(subnacontigsshort_sam, True, nbins)
    subnacontigslong_nalbyerr = extract_from_bam(subnacontigslong_sam, False, nbins)
    subnacontigs_nalbyerr = subnacontigsshort_nalbyerr + subnacontigslong_nalbyerr
    print("On NA contigs (sub):", len(subnacontigs_nalbyerr))

    # Plotting
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey = True, tight_layout=True, figsize=(13, 4))

    ax1.axhline(tot_strings, linewidth=1, color='r', label="Total specific")
    ax1.hist(hgcontigs_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = True,
             align = "left", histtype="step",
             color="seagreen", alpha=0.75,
             label="Cumulative count")
    ax1.hist(hgcontigs_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = False,
             align = "left", histtype="stepfilled",
             color="seagreen", alpha=0.75,
             label="Count")
    ax1.set_xticks(np.arange(0, nbins))
    ax1.set_title("(a) HG00733 contigs")
    ax1.legend(loc=4, fancybox=True, shadow=True, fontsize="small", bbox_to_anchor=(1, 0.07))

    ax2.axhline(tot_strings, linewidth=1, color='r', label="Total specific")
    ax2.hist(nacontigs_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = True,
             align = "left", histtype="step",
             color="seagreen", alpha=0.75,
             label="Cumulative count")
    ax2.hist(nacontigs_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = False,
             align = "left", histtype="stepfilled",
             color="seagreen", alpha=0.75,
             label="Count")
    ax2.set_xticks(np.arange(0, nbins))
    ax2.set_title("(b) NA19240 contigs")
    ax2.legend(loc=4, fancybox=True, shadow=True, fontsize="small", bbox_to_anchor=(1, 0.07))

    ax3.hist([haplocomp, nonhaplocomp],
             bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
             histtype="bar", cumulative=False,
             color=["steelblue", "darkorchid"], label=["Haplo-comp.", "Non haplo-comp."])
    ax3.set_title("(c) HG00733 haplotypes")
    ax3.set_xticks(np.arange(0, nbins))
    ax3.legend(loc=1, fancybox=True, shadow=True, fontsize="small")

    ax4.hist([subhgcontigs_nalbyerr, subnacontigs_nalbyerr],
             bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
             histtype="bar", cumulative=False,
             color=["plum", "pink"], label=["HG00733", "NA19240"])
    ax4.set_title("(d) Contigs (non haplo-comp.)")
    ax4.set_xticks(np.arange(0, nbins))
    ax4.legend(loc=1, fancybox=True, shadow=True, fontsize="small")

    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
    plt.xticks()

    plt.xlabel("# Base Differences")
    plt.ylabel("# Primary Alignments") #, labelpad = 17)

    plt.savefig(out_path, dpi=300)
    # plt.show()

if __name__ == "__main__":
    main()
