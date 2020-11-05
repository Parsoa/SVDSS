import sys

from Bio import SeqIO
import pysam
from pysam import VariantFile

###########################
### ALLELE-BASED RECALL ###
###########################
def get_atype(alen, splitL):
    if alen == 1:
        return "SNP"
    elif alen < splitL:
        return "<"
    else:
        return ">="

def extract_alleles(bed_path, splitL):
    alleles = {}
    for line in open(bed_path, 'r'):
        line = line.strip('\n').split('\t')
        aidx, alen, gt, covering = line[3], int(len(line[5])), int(line[-2]), int(line[-1])
        atype = get_atype(alen, splitL)
        hap = '1' if aidx.endswith('1') else '2'
        gt = "Ref" if gt == 0 else "Alt"
        alleles[aidx] = (covering>0, atype, hap, gt)
    return alleles

def recall():
    # TODO this can be done better (more pythonic)
    bed1_path = sys.argv[1]    # bedintersect -c between alignments and unique variants on first haplotype
    bed2_path = sys.argv[2]    # bedintersect -c between alignments and unique variants on second haplotype
    splitL = int(sys.argv[3])  # split length (indels/SVs)

    alleles1 = extract_alleles(bed1_path, splitL)
    alleles2 = extract_alleles(bed2_path, splitL)

    alleles = {**alleles1, **alleles2}

    found = {"SNP" : {"Ref" : 0, "Alt" : 0}, "<" : {"Ref" : 0, "Alt" : 0}, ">=" : {"Ref" : 0, "Alt" : 0}}
    tot = {"SNP" : {"Ref" : 0, "Alt" : 0}, "<" : {"Ref" : 0, "Alt" : 0}, ">=" : {"Ref" : 0, "Alt" : 0}}
    for aidx, (covering, atype, hap, gt) in alleles.items():
        tot[atype][gt] += 1
        if covering:
            found[atype][gt] += 1

    print(f"Split at", splitL)
    print("")
    print("Type", "GT", "Found", "Total", "%", sep='\t')
    print("---")
    ntot = {"Ref": 0, "Alt": 0}
    nfound = {"Ref": 0, "Alt": 0}
    for idx in ["SNP", "<", ">="]:
        for gt in ["Ref", "Alt"]:
            print(idx, gt, found[idx][gt], tot[idx][gt], round(found[idx][gt]/tot[idx][gt]*100, 2), sep='\t')
            nfound[gt] += found[idx][gt]
            ntot[gt] += tot[idx][gt]
    print("---")
    print("All", "Ref", nfound["Ref"], ntot["Ref"], round(nfound["Ref"]/ntot["Ref"]*100, 2), sep='\t')
    print("All", "Alt", nfound["Alt"], ntot["Alt"], round(nfound["Alt"]/ntot["Alt"]*100, 2), sep='\t')



##########################
# CONTIG-BASED PRECISION #
##########################
def contigbased_precision():
    bampath1 = sys.argv[1]
    bampath2 = sys.argv[2]
    fq_path = sys.argv[3]

    total_specific = 0
    for record in SeqIO.parse(fq_path, "fastq"):
        total_specific += 1

    perfect1 = set()
    bamfile1 = pysam.AlignmentFile(bampath1, "rb")
    for al in bamfile1.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        if bad_bases + deletions == 0:
            perfect1.add(al.query_name)

    perfect2 = set()
    bamfile2 = pysam.AlignmentFile(bampath2, "rb")
    for al in bamfile2.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        if bad_bases + deletions == 0:
            perfect2.add(al.query_name)

    print("Total (T):", total_specific, sep='\t')
    print("Perfect1 (P1):", len(perfect1), sep='\t')
    print("Perfect2 (P2):", len(perfect2), sep='\t')
    print("P1 & P2:", len(perfect1 & perfect2), sep='\t')
    print("P1 - P2:", len(perfect1 - perfect2), sep='\t')
    print("(P1-P2)/T:", round(len(perfect1 - perfect2)/total_specific*100, 2), sep='\t')



################################
### HAPLOTYPE-BASE PRECISION ###
################################
def parse_waobed(bed_path, firsthap = True):
    alignments = {}
    all_alignments = {}
    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        ref, start, end, ridx, vidx, covering = line[0], line[1], line[2], line[3], line[9], line[-1]
        if ref.endswith("_2") and firsthap:
            continue
        posidx = ridx + ":" + start + "-" + end
        if posidx not in alignments:
            alignments[posidx] = set()
        if covering != '0':
            alignments[posidx].add(int(vidx.split('_')[0]))
    return alignments

def parse_vcf(vcf_path, sample1 = "HG00733", sample2 = "NA19240"):
    variants = {}
    vcf = VariantFile(vcf_path)
    idx = 1
    for record in vcf.fetch():
        gt1, gt2 = record.samples[sample1]["GT"]
        gt1_2, gt2_2 = record.samples[sample2]["GT"]
        variants[idx] = [(gt1, gt2), (gt1_2, gt2_2)] # maybe we can use 4 lists
        idx += 1
    return variants

def parse_outsiders(bed_path):
    outsiders_dict = {}
    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        ridx, flag = line[3], line[-1]
        if ridx not in outsiders_dict:
            outsiders_dict[ridx] = False
        outsiders_dict[ridx] = outsiders_dict[ridx] or flag == '0'
    outsiders = set()
    for ridx,flag in outsiders_dict.items():
        if not flag:
            outsiders.add(ridx)
    print(len(outsiders), len(outsiders_dict), len(outsiders)/len(outsiders_dict) if len(outsiders_dict) > 0 else 0)
    return outsiders

def check_haplo_uniqueness(variants, firsthap = True):
    HGhaplo = [v[0][not firsthap] for v in variants]
    NAhaplo1 = [v[1][0] for v in variants]
    NAhaplo2 = [v[1][1] for v in variants]
    return HGhaplo != NAhaplo1 and HGhaplo != NAhaplo2

def precision():
    vcf_path = sys.argv[1]  # all variants
    bed1_path = sys.argv[2] # bedintersect -wao between alignments and variants on first haplotype
    bed2_path = sys.argv[3] # bedintersect -wao between alignments and variants on second haplotype
    bed3_path = sys.argv[4] # bedintersect -c between alignments and regions not called by HGSV consortium

    variants = parse_vcf(vcf_path)

    alignments1 = parse_waobed(bed1_path, True)
    alignments2 = parse_waobed(bed2_path, False)

    outsiders = parse_outsiders(bed3_path)

    covering_results = {}
    for posidx,vidxs in alignments1.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, True)
            covering_results[ridx] = covering_results[ridx] or unique

    for posidx,vidxs in alignments2.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, False)
            covering_results[ridx] = covering_results[ridx] or unique

    total = 0
    covering = 0
    outside = 0
    for ridx, flag in covering_results.items():
        total += 1
        if flag:
            covering += 1
        else:
            if ridx in outsiders:
                outside += 1
    print("Covering:", covering)
    print("Outsiders:", outside)
    print("Total:", total)
    print("P:", round(covering/total*100, 3))
    print("P (w/o outsiders):", round(covering/(total-outside)*100, 3))

# I had to duplicate the function to avoid a cycle in snakemake
def extract_uncovering():
    vcf_path = sys.argv[1]  # all variants
    bed1_path = sys.argv[2] # bedintersect -wao between alignments and variants on first haplotype
    bed2_path = sys.argv[3] # bedintersect -wao between alignments and variants on second haplotype

    variants = parse_vcf(vcf_path)

    alignments1 = parse_waobed(bed1_path, True)
    alignments2 = parse_waobed(bed2_path, False)

    covering_results = {}
    for posidx,vidxs in alignments1.items():
        ridx = posidx.split(':')[0]
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, True)
            covering_results[ridx] = covering_results[ridx] or unique

    for posidx,vidxs in alignments2.items():
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

if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "pre":
        precision()
    elif mode == "pre1":
        contigbased_precision()
    elif mode == "rec":
        recall()
    elif mode == "exunc":
        extract_uncovering()
