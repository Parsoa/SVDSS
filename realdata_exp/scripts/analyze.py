import sys

from Bio import SeqIO
import pysam

def extract_reclist():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        s, e, idx, l, count = int(line[1]), int(line[2]), line[3], int(line[6]), int(line[-1])
        l1 = e - s - 1
        if l < 0:
            l1 = -l1
        l = l1
        variants[idx] = (count, -1, l)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        s, e, idx, l, count = int(line[1]), int(line[2]), line[3], int(line[6]), int(line[-1])
        l1 = e - s - 1
        if l < 0:
            l1 = -l1
        l = l1
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

def recall():
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

def new_precision():
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
        new_precision()
    elif mode == "reclist":
        extract_reclist()
    elif mode == "rec":
        recall()
    elif mode == "exunc":
        extract_uncovering()
