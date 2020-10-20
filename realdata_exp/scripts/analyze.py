import sys

from Bio import SeqIO
import pysam

def get_atype(alen, l):
    if alen == 1:
        return "SNP"
    elif alen < l:
        return "<"
    else:
        return ">="

def extract_alleles(bed_path, alt):
    alleles = {}
    for line in open(bed_path, 'r'):
        line = line.strip('\n').split('\t')
        aidx, alen, gt, covering = line[3], int(len(line[5])), int(line[-2]), int(line[-1])
        atype = get_atype(alen, alt)
        hap = '1' if aidx.endswith('1') else '2'
        gt = "Ref" if gt == 0 else "Alt"
        alleles[aidx] = (covering>0, atype, hap, gt)
    return alleles

def recall():
    bed1_path = sys.argv[1]
    bed2_path = sys.argv[2]
    alt = int(sys.argv[3])

    alleles1 = extract_alleles(bed1_path, alt)
    alleles2 = extract_alleles(bed2_path, alt)

    alleles = {**alleles1, **alleles2}

    found = {"SNP" : {"Ref" : 0, "Alt" : 0}, "<" : {"Ref" : 0, "Alt" : 0}, ">=" : {"Ref" : 0, "Alt" : 0}}
    tot = {"SNP" : {"Ref" : 0, "Alt" : 0}, "<" : {"Ref" : 0, "Alt" : 0}, ">=" : {"Ref" : 0, "Alt" : 0}}
    for aidx, (covering, atype, hap, gt) in alleles.items():
        tot[atype][gt] += 1
        if covering:
            found[atype][gt] += 1

    print(f"Split at", alt)
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

def extract_reclist():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        s, e, idx, l, count = int(line[1]), int(line[2]), line[3], int(line[6]), int(line[-1])
        variants[idx] = (count, -1, l)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        s, e, idx, l, count = int(line[1]), int(line[2]), line[3], int(line[6]), int(line[-1])
        if idx in variants:
            variants[idx] = (variants[idx][0], count, l)
        else:
            variants[idx] = (-1, count, l)

    for idx, (c1,c2,l) in variants.items():
        ov = 0
        if c1 <= 0 and c2 <= 0:
            ov = 0
        elif c1 > 0 and c2 > 0:
            ov = 3
        elif c1 > 0 and c2 <= 0:
            ov = 1
        elif c1 <= 0 and c2 > 0:
            ov = 2
        print(idx, l, ov)

def recall_from_reclist():
    fpath = sys.argv[1]
    l_thresh = int(sys.argv[2])

    found = {"SNP" : 0, "<" : 0, ">=" : 0}
    tot = {"SNP" : 0, "<" : 0, ">=" : 0}
    for line in open(fpath):
        idx, l, ov = line.strip('\n').split(' ')
        l, ov = abs(int(l)), int(ov)

        idx = "SNP"
        if l == 0:
            idx = "SNP"
        elif l < l_thresh:
            idx = "<"
        else:
            idx = ">="

        tot[idx] += 1
        if ov > 0:
            found[idx] += 1

    print(f"Split at", l_thresh)
    print("")
    print("Type", "Found", "Total", "%", sep='\t')
    print("---")
    for idx in ["SNP", "<", ">="]:
        print(idx, found[idx], tot[idx], round(found[idx]/tot[idx]*100, 2), sep='\t')
    print("---")
    print("All", sum(found.values()), sum(tot.values()), round(sum(found.values())/sum(tot.values())*100, 2), sep='\t')

def precision():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if idx not in variants:
            variants[idx] = (-1,-1)
        variants[idx] = (max(count, variants[idx][0]), -1)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if idx not in variants:
            variants[idx] = (-1,-1)
        variants[idx] = (variants[idx][0], max(count, variants[idx][1]))

    none = 0
    both = 0
    only1 = 0
    only2 = 0
    total = 0
    for (c1,c2) in variants.values():
        total += 1
        if c1 <= 0 and c2 <= 0:
            none += 1
        elif c1 > 0 and c2 > 0:
            both += 1
        elif c1 > 0 and c2 <= 0:
            only1 += 1
        elif c1 <= 0 and c2 > 0:
            only2 += 1

    print("None", "Both", "1", "2", "Total", "Precision", sep=',')
    print(none, both, only1, only2, total, round(((both+only1+only2)/total)*100, 2), sep=',')

def extract_uncovering():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]
    fq_path = sys.argv[3] if len(sys.argv) == 4 else ""

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if idx not in variants:
            variants[idx] = (-1,-1)
        variants[idx] = (max(count, variants[idx][0]), -1)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if idx not in variants:
            variants[idx] = (-1,-1)
        variants[idx] = (variants[idx][0], max(count, variants[idx][1]))

    if fq_path == "":
        for idx, (c1,c2) in variants.items():
            if c1 <= 0 and c2 <= 0:
                print(idx)
    else:
        # FASTQ OUTPUT
        uncovering_idxs = set()
        for idx, (c1,c2) in variants.items():
            if c1 <= 0 and c2 <= 0:
                uncovering_idxs.add(idx)
        for record in SeqIO.parse(fq_path, "fastq"):
            if record.id in uncovering_idxs:
                SeqIO.write(record, sys.stdout, "fastq")

def new_contigbased_precision():
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

if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "pre":
        precision()
    elif mode == "newpre":
        new_contigbased_precision()
    elif mode == "reclist":
        extract_reclist()
    elif mode == "rec":
        recall()
    elif mode == "recfromlist":
        recall_from_reclist()
    elif mode == "exunc":
        extract_uncovering()
