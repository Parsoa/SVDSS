import sys

from Bio import SeqIO
import pysam
from pysam import VariantFile

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

###########################
### ALLELE-BASED RECALL ###
###########################


def get_vtype(idx, splitL):
    idx = idx.split('-')
    if idx[2] == "SNV":
        return "SNP"
    else:
        if int(idx[-1]) < splitL:
            return "INDEL"
        else:
            return "SV"


def extract_alleles(bed_path, splitL):
    alleles = {}
    for line in open(bed_path, 'r'):
        line = line.strip('\n').split('\t')
        chrom, idx, gt, covering = line[0], line[3], int(
            line[-2]), int(line[-1])
        vtype = get_vtype(idx, splitL)
        hap = '2' if chrom.endswith("_2") else '1'
        idx += hap
        gt = "Ref" if gt == 0 else "Alt"
        if idx in alleles:
            alleles[idx] = (covering > 0 or alleles[idx][0], vtype, hap, gt)
        else:
            alleles[idx] = (covering > 0, vtype, hap, gt)
    return alleles


def recall():
    # TODO this can be done better (more pythonic)
    # bedintersect -c between alignments and unique variants on first haplotype
    bed1_path = sys.argv[1]
    # bedintersect -c between alignments and unique variants on second haplotype
    bed2_path = sys.argv[2]
    splitL = int(sys.argv[3])

    alleles1 = extract_alleles(bed1_path, splitL)
    alleles2 = extract_alleles(bed2_path, splitL)

    alleles = {**alleles1, **alleles2}

    found = {"SNP": {"Ref": 0, "Alt": 0}, "INDEL": {
        "Ref": 0, "Alt": 0}, "SV": {"Ref": 0, "Alt": 0}}
    tot = {"SNP": {"Ref": 0, "Alt": 0}, "INDEL": {
        "Ref": 0, "Alt": 0}, "SV": {"Ref": 0, "Alt": 0}}
    for aidx, (covering, vtype, hap, gt) in alleles.items():
        tot[vtype][gt] += 1
        if covering:
            found[vtype][gt] += 1

    print(f"Split at", splitL)
    print("")
    print("Type", "GT", "Found", "Total", "%", sep='\t')
    print("---")
    ntot = {"Ref": 0, "Alt": 0}
    nfound = {"Ref": 0, "Alt": 0}
    for idx in ["SNP", "INDEL", "SV"]:
        for gt in ["Ref", "Alt"]:
            print(idx, gt, found[idx][gt], tot[idx][gt], round(
                found[idx][gt]/tot[idx][gt]*100, 2) if tot[idx][gt] != 0 else 0, sep='\t')
            nfound[gt] += found[idx][gt]
            ntot[gt] += tot[idx][gt]
    print("---")
    print("All", "Ref", nfound["Ref"], ntot["Ref"], round(
        nfound["Ref"]/ntot["Ref"]*100 if ntot["Ref"] != 0 else 0, 2), sep='\t')
    print("All", "Alt", nfound["Alt"], ntot["Alt"], round(
        nfound["Alt"]/ntot["Alt"]*100, 2) if ntot["Alt"] else 0, sep='\t')

#################################
### HAPLOTYPE-BASED PRECISION ###
#################################


def count_sample(fq_path):
    n = 0
    den = 0
    for line in open(fq_path):
        if den == 0:
            den = 2 if line.startswith('>') else 4
        n += 1
    return n/den if den > 0 else 0


def parse_waobed(bed_path, firsthap=True):
    alignments = {}
    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        ref, start, end, ridx, vidx, covering = line[0], line[1], line[2], line[3], line[9], line[-1]
        if ref.endswith("_2") and firsthap:
            continue
        posidx = ridx + ":" + start + "-" + end # we need positions to cluster close variants (we will "reconstruct" the haplotype to see if the portion covered by the alignment is specific)
        if posidx not in alignments:
            alignments[posidx] = set()
        if covering != '0':
            alignments[posidx].add(vidx)
    return alignments


def parse_vcf(vcf_path, sample1="HG00733", sample2="NA19240"):
    variants = {}
    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        idx = record.id
        gt1, gt2 = record.samples[sample1]["GT"]
        gt1_2, gt2_2 = record.samples[sample2]["GT"]
        # maybe we can use 4 lists
        variants[idx] = [(gt1, gt2), (gt1_2, gt2_2)]
    return variants


def check_haplo_uniqueness(variants, firsthap=True):
    HGhaplo = [v[0][not firsthap] for v in variants]
    NAhaplo1 = [v[1][0] for v in variants]
    NAhaplo2 = [v[1][1] for v in variants]
    return HGhaplo != NAhaplo1 and HGhaplo != NAhaplo2


def precision():
    fq_path = sys.argv[1]   # sample
    vcf_path = sys.argv[2]  # all variants
    # bedintersect -wao between alignments and variants on first haplotype
    bed1_path = sys.argv[3]
    # bedintersect -wao between alignments and variants on second haplotype
    bed2_path = sys.argv[4]

    total = count_sample(fq_path)

    variants = parse_vcf(vcf_path)

    alignments1 = parse_waobed(bed1_path, True)
    alignments2 = parse_waobed(bed2_path, False)

    covering_results = {}
    for posidx, vidxs in alignments1.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, True)
            covering_results[ridx] = covering_results[ridx] or unique

    for posidx, vidxs in alignments2.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, False)
            covering_results[ridx] = covering_results[ridx] or unique

    covering = 0
    for ridx, flag in covering_results.items():
        if flag:
            covering += 1

    print("Covering:", covering)
    print("Total:", total)
    print("P:", round(covering/total*100 if total != 0 else 0, 3))


def extract_uncovered():
    bedpath1 = sys.argv[1]
    bedpath2 = sys.argv[2]
    vcfpath = sys.argv[3]
    outprefix = sys.argv[4]
    splitL = 250

    alleles1 = extract_alleles(bedpath1, splitL)
    alleles2 = extract_alleles(bedpath2, splitL)

    alleles = {**alleles1, **alleles2}
    fns = set()
    for aidx, (covering, atype, hap, gt) in alleles.items():
        if not covering:
            fn = int(aidx.split('_')[0])
            fns.add(fn)

    vcf = VariantFile(vcfpath)
    smallout = open(outprefix + ".short.vcf", 'w')
    midout = open(outprefix + ".mid.vcf", 'w')
    longout = open(outprefix + ".long.vcf", 'w')
    print('\n'.join(str(vcf.header).split('\n')), end='', file=smallout)
    print('\n'.join(str(vcf.header).split('\n')), end='', file=midout)
    print('\n'.join(str(vcf.header).split('\n')), end='', file=longout)
    n = 1
    for record in vcf.fetch():
        if n in fns:
            idx = record.id.split('-')
            if idx[2] == "SNV":
                print(record, end='', file=smallout)
            else:
                if int(idx[-1]) < splitL:
                    print(record, end='', file=midout)
                else:
                    print(record, end='', file=longout)
        n += 1
    smallout.close()
    midout.close()
    longout.close()

##########################
# CONTIG-BASED PRECISION #
##########################


def contigbased_precision():
    if len(sys.argv) == 4:
        bampath1 = sys.argv[1]  # all on hg
        bampath2 = ""
        bampath3 = sys.argv[2]  # all on na
        bampath4 = ""
        fq_path = sys.argv[3]   # sample
    elif len(sys.argv) == 6:
        bampath1 = sys.argv[1]  # short on hg
        bampath2 = sys.argv[2]  # long on hg
        bampath3 = sys.argv[3]  # short on na
        bampath4 = sys.argv[4]  # long on na
        fq_path = sys.argv[5]   # sample

    total_specific = count_sample(fq_path)

    # counting perfect on hg
    perfect1 = set()
    bamfile1 = pysam.AlignmentFile(bampath1, "rb")  # bbmap
    for al in bamfile1.fetch():
        # good_bases = al.get_cigar_stats()[0][7]


        if al.is_supplementary or al.is_unmapped or al.is_secondary:
            continue
        good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")


        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        if bad_bases + deletions == 0:
            perfect1.add(al.query_name)

    if bampath2 != "":
        bamfile2 = pysam.AlignmentFile(bampath2, "rb")  # minimap2
        for al in bamfile2.fetch():
            if al.is_supplementary or al.is_unmapped or al.is_secondary:
                continue
            good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
            bad_bases = al.query_length - good_bases
            deletions = al.get_cigar_stats()[0][2]
            if bad_bases + deletions == 0:
                perfect1.add(al.query_name)

    # counting perfect on na
    perfect2 = set()
    bamfile3 = pysam.AlignmentFile(bampath3, "rb")  # bbmap
    for al in bamfile3.fetch():
        # good_bases = al.get_cigar_stats()[0][7]


        if al.is_supplementary or al.is_unmapped or al.is_secondary:
            continue
        good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")


        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        if bad_bases + deletions == 0:
            perfect2.add(al.query_name)

    if bampath4 != "":
        bamfile4 = pysam.AlignmentFile(bampath4, "rb")  # minimap2
        for al in bamfile4.fetch():
            if al.is_supplementary or al.is_unmapped or al.is_secondary:
                continue
            good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
            bad_bases = al.query_length - good_bases
            deletions = al.get_cigar_stats()[0][2]
            if bad_bases + deletions == 0:
                perfect2.add(al.query_name)

    print("Total (T):", total_specific, sep='\t')
    print("Perfect1 (P1):", len(perfect1), sep='\t')
    print("Perfect2 (P2):", len(perfect2), sep='\t')
    print("P1 & P2:", len(perfect1 & perfect2), sep='\t')
    print("P1 - P2:", len(perfect1 - perfect2), sep='\t')
    print("(P1-P2)/T:", round(len(perfect1 - perfect2) /
                              total_specific*100, 2), sep='\t')


def cplot():
    ### quality of alignments to contigs ###
    bampath1 = sys.argv[1]  # short on hg
    bampath2 = sys.argv[2]  # long on hg
    bampath3 = sys.argv[3]  # short on na
    bampath4 = sys.argv[4]  # long on na
    fq_path = sys.argv[5]  # sample
    out_path = sys.argv[6]  # sample

    tot_strings = count_sample(fq_path)
    print("Tot strings: ", tot_strings)

    # counting perfect on hg
    data1 = []
    bamfile1 = pysam.AlignmentFile(bampath1, "rb")  # bbmap
    for al in bamfile1.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data1.append(errors)

    bamfile2 = pysam.AlignmentFile(bampath2, "rb")  # minimap2
    for al in bamfile2.fetch():
        if al.is_supplementary or al.is_unmapped or al.is_secondary:
            continue
        good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data1.append(errors)
    print("bam1 limits: ", min(data1), max(data1))

    # counting perfect on na
    data2 = []
    bamfile3 = pysam.AlignmentFile(bampath3, "rb")  # bbmap
    for al in bamfile3.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data2.append(errors)

    bamfile4 = pysam.AlignmentFile(bampath4, "rb")  # minimap2
    for al in bamfile4.fetch():
        if al.is_supplementary or al.is_unmapped or al.is_secondary:
            continue
        good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions
        data2.append(errors)

    print("bam2 limits: ", min(data2), max(data2))
    print(len(data1), len(data2), len(data1)+len(data2))

    nbins = 7
    fig, (ax1, ax2) = plt.subplots(
        1, 2, sharey=True, tight_layout=True, figsize=(13, 7))

    ax1.axhline(tot_strings, 0, 100, linewidth=1,
                color='r', label="Total specific")
    ax1.hist(data1, bins=np.arange(0, nbins+1), range=(0, nbins+1), align="left",
             color="seagreen", alpha=0.75, cumulative=True, histtype="step", label="Cumulative count")
    ax1.hist(data1, bins=np.arange(0, nbins+1), range=(0, nbins+1),
             align="left", color="seagreen", histtype="stepfilled", label="Count")
    ax1.set_title("HG00733")
    ax1.set_ylabel("# Primary Alignments")

    ax2.axhline(tot_strings, 0, 100, linewidth=1,
                color='r', label="Total specific")
    ax2.hist(data2, bins=np.arange(0, nbins+1), range=(0, nbins+1), align="left",
             color="seagreen", alpha=0.75, cumulative=True, histtype="step", label="Cumulative count")
    ax2.hist(data2, bins=np.arange(0, nbins+1), range=(0, nbins+1),
             align="left", color="seagreen", histtype="stepfilled", label="Count")
    ax2.set_title("NA19240")

    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor='none', top=False,
                    bottom=False, left=False, right=False)
    plt.xticks()
    plt.xlabel("# Base Differences")

    ax2.legend(loc=4, fancybox=True, shadow=True)
    plt.savefig(out_path)
    # plt.show()

# I had to duplicate the function to avoid a cycle in snakemake


def extract_nonhaplo():
    vcf_path = sys.argv[1]  # all variants
    # bedintersect -wao between alignments and variants on first haplotype
    bed1_path = sys.argv[2]
    # bedintersect -wao between alignments and variants on second haplotype
    bed2_path = sys.argv[3]

    variants = parse_vcf(vcf_path)

    alignments1 = parse_waobed(bed1_path, True)
    alignments2 = parse_waobed(bed2_path, False)

    covering_results = {}
    for posidx, vidxs in alignments1.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, True)
            covering_results[ridx] = covering_results[ridx] or unique

    for posidx, vidxs in alignments2.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, False)
            covering_results[ridx] = covering_results[ridx] or unique

    for ridx, flag in covering_results.items():
        if not flag:
            print(ridx)


def plot_nonhaplo():
    uncovering_path = sys.argv[1]
    haplosam_path = sys.argv[2]
    contigssam1_path = sys.argv[3]
    contigssam2_path = sys.argv[4]
    out_path = sys.argv[5]

    # plot (hardcoded) parameters
    nbins = 7
    barw = 0.35
    Xs = list(range(0, nbins))

    # Parsing uncoverings
    uncovering_idxs = set()
    for line in open(uncovering_path):
        uncovering_idxs.add(line.strip('\n'))

    # Parsing haplo sam (left plot)
    covering = []
    not_covering = []
    haplosam = pysam.AlignmentFile(haplosam_path, 'r')
    for al in haplosam.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        # bbmap:
        good_bases = al.get_cigar_stats()[0][7]
        # minimap2:
        # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        errors = bad_bases + deletions

        if al.query_name in uncovering_idxs:
            not_covering.append(errors)
        else:
            covering.append(errors)

    print("Covering:", len(covering))
    print("Limits:", min(covering), max(covering))
    print("---")
    print("Not covering:", len(not_covering))
    print("Limits:", min(not_covering), max(not_covering))

    # Parsing contigs sams (right plot)
    data1 = []
    contigssam1 = pysam.AlignmentFile(contigssam1_path, 'r')
    for al in contigssam1.fetch():
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
    contigssam2 = pysam.AlignmentFile(contigssam2_path, 'r')
    for al in contigssam2.fetch():
        if al.is_secondary or al.is_unmapped or al.is_supplementary:
            continue
        # bbmap:
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
    fig, (ax1, ax2) = plt.subplots(
        1, 2, sharey=True, tight_layout=True, figsize=(13, 7))

    ax1.hist([covering, not_covering],
             bins=np.arange(0, nbins+1), range=(0, nbins+1), align="left",
             histtype="bar", cumulative=False,
             color=["steelblue", "darkorchid"], label=["Covering", "Not covering"])
    ax1.set_ylabel("# Primary Alignments")
    ax1.legend()

    ax2.bar(Xs, [data1_dict[x] if x in data1_dict else 0 for x in Xs],
            color="plum", width=-barw, align="edge", label="HG00733")
    ax2.bar(Xs, [data2_dict[x] if x in data2_dict else 0 for x in Xs],
            color="pink", width=barw, align="edge", label="NA19240")
    ax2.legend()

    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", top=False,
                    bottom=False, left=False, right=False)
    plt.xticks()
    plt.xlabel("# Base Differences")

    plt.savefig(out_path)

### DATA ANALYSIS PLOTS ###
def plot_samplelengths():
    fq_path1 = sys.argv[1]
    fq_path2 = sys.argv[2]
    out_path = sys.argv[3]

    data1 = []
    with open(fq_path1, 'r') as fq:
        for i,line in enumerate(fq):
            if i-1 % 4 == 0:
                l = len(line)-1
                data1.append(l)
    data2 = []
    with open(fq_path2, 'r') as fq:
        for i,line in enumerate(fq):
            if i-1 % 4 == 0:
                l = len(line)-1
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

if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "pre":
        precision()
    elif mode == "pre1":
        contigbased_precision()
    elif mode == "rec":
        recall()
    elif mode == "fn":
        extract_uncovered()
    elif mode == "nonhaplo":
        extract_nonhaplo()
    elif mode == "plotnonhaplo":
        plot_nonhaplo()
    elif mode == "cplot":
        cplot()
    elif mode == "sampleshist":
        plot_samplelengths()
    elif mode == "specorr":
        plot_correlation()
