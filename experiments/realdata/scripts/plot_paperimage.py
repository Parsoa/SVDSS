import sys

from Bio import SeqIO
import pysam

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_al_errors(al):
    # bbmap:
    good_bases = al.get_cigar_stats()[0][7]
    # minimap2:
    # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
    bad_bases = al.query_length - good_bases
    deletions = al.get_cigar_stats()[0][2]
    return bad_bases + deletions

def extract_from_bam(bampath, cutoff = float("inf")):
    data = []
    bamfile = pysam.AlignmentFile(bampath, "rb")
    for al in bamfile.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        errors = get_al_errors(al)
        if errors < cutoff:
            data.append(errors)
    return data

def main():
    nbins = 6
    
    fq_path = sys.argv[1]
    HGcontigssam_path = sys.argv[2]
    NAcontigssam_path = sys.argv[3]
    nonhaplo_path = sys.argv[4]
    HGhaplosam_path = sys.argv[5]
    HGcontigssubsam1_path = sys.argv[6]
    NAcontigssubsam1_path = sys.argv[7]
    out_path = sys.argv[8]

    tot_strings = 0
    # for record in SeqIO.parse(fq_path, "fastq"):
    #     tot_strings += 1
    tot_strings = 7041081
    print("Tot strings: ", tot_strings)

    HGcon_nalbyerr = extract_from_bam(HGcontigssam_path, nbins)
    print("HGcontigs:", len(HGcon_nalbyerr))
    NAcon_nalbyerr = extract_from_bam(NAcontigssam_path, nbins)
    print("NAcontigs:", len(NAcon_nalbyerr))

    # Parsing non-haplo-aware strings (after alignment to HG haplo)
    nonhaplo_idxs = set()
    for line in open(nonhaplo_path):
        nonhaplo_idxs.add(line.strip('\n'))
    print("Non Haplo-Aware:", len(nonhaplo_idxs))

    # Parsing haplo sam
    haploaware = []
    nonhaploaware = []
    haplosam = pysam.AlignmentFile(HGhaplosam_path, 'r')
    for al in haplosam.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        errors = get_al_errors(al)
        if errors < nbins:
            if al.query_name in nonhaplo_idxs:
                nonhaploaware.append(errors)
            else:
                haploaware.append(errors)
    print("Haplo-Aware:", len(haploaware))
    print("Non Haplo-Aware:", len(nonhaploaware))

    HGconsub_nalbyerr = extract_from_bam(HGcontigssubsam1_path, nbins)
    print("HGcontigs (sub):", len(HGconsub_nalbyerr))
    NAconsub_nalbyerr = extract_from_bam(NAcontigssubsam1_path, nbins)
    print("NAcontigs (sub):", len(NAconsub_nalbyerr))

    # Plotting
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey = True, tight_layout=True, figsize=(13, 4))

    ax1.axhline(tot_strings, linewidth=1, color='r', label="Total specific")    
    ax1.hist(HGcon_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = True,
             align = "left", histtype="step",
             color="seagreen", alpha=0.75,
             label="Cumulative count")
    ax1.hist(HGcon_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = False,
             align = "left", histtype="stepfilled",
             color="seagreen", alpha=0.75,
             label="Count")
    ax1.set_xticks(np.arange(0, nbins))
    ax1.set_title("(a) HG00733 contigs")
    ax1.legend(loc=4, fancybox=True, shadow=True, fontsize="small", bbox_to_anchor=(1, 0.07))
    
    ax2.axhline(tot_strings, linewidth=1, color='r', label="Total specific")
    ax2.hist(NAcon_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = True,
             align = "left", histtype="step",
             color="seagreen", alpha=0.75,
             label="Cumulative count")
    ax2.hist(NAcon_nalbyerr,
             bins=np.arange(0, nbins+1), range=(0,nbins+1),
             cumulative = False,
             align = "left", histtype="stepfilled",
             color="seagreen", alpha=0.75,
             label="Count")
    ax2.set_xticks(np.arange(0, nbins))
    ax2.set_title("(b) NA19240 contigs")
    ax2.legend(loc=4, fancybox=True, shadow=True, fontsize="small", bbox_to_anchor=(1, 0.07))

    ax3.hist([haploaware, nonhaploaware],
             bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
             histtype="bar", cumulative=False,
             color=["steelblue", "darkorchid"], label=["Haplo-aware", "Non haplo-aware"])
    ax3.set_title("(c) HG00733 haplotypes")
    ax3.set_xticks(np.arange(0, nbins))
    ax3.legend(loc=1, fancybox=True, shadow=True, fontsize="small")

    ax4.hist([HGconsub_nalbyerr, NAconsub_nalbyerr],
             bins=np.arange(0, nbins+1), range=(0,nbins+1), align = "left",
             histtype="bar", cumulative=False,
             color=["plum", "pink"], label=["HG00733", "NA19240"])
    ax4.set_title("(d) Contigs (non haplo-aware)")
    ax4.set_xticks(np.arange(0, nbins))
    ax4.legend(loc=1, fancybox=True, shadow=True, fontsize="small")

    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
    plt.xticks()
    
    plt.xlabel("# Base Differences")
    plt.ylabel("# Primary Alignments") #, labelpad = 17)

    plt.savefig(out_path)
    # plt.show()
    
if __name__ == "__main__":
    main()
