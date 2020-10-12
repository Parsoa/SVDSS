import sys

from Bio import SeqIO

def recall():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, l, count = line[3], int(line[6]), int(line[-1])
        variants[idx] = (count, -1, l)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, l, count = line[3], int(line[6]), int(line[-1])
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

def precision():
    bed_path_1 = sys.argv[1]
    bed_path_2 = sys.argv[2]

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        variants[idx] = (count, -1)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if idx in variants:
            variants[idx] = (variants[idx][0], count)
        else:
            variants[idx] = (-1, count)

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
    # fq_path = sys.argv[3]

    variants = {}
    for line in open(bed_path_1):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        variants[idx] = (count, -1)

    for line in open(bed_path_2):
        line = line.strip('\n').split('\t')
        idx, count = line[3], int(line[-1])
        if idx in variants:
            variants[idx] = (variants[idx][0], count)
        else:
            variants[idx] = (-1, count)

    for idx, (c1,c2) in variants.items():
        if c1 <= 0 and c2 <= 0:
            print(idx)

    ### IF FASTQ OUTPUT
    # uncovering_idxs = set()
    # for idx, (c1,c2) in variants.items():
    #     if c1 <= 0 and c2 <= 0:
    #         uncovering_idxs.add(idx)
    # for record in SeqIO.parse(fq_path, "fastq"):
    #     if record.id in uncovering_idxs:
    #         SeqIO.write(record, sys.stdout, "fastq")

def new_prec():
    bampath1 = sys.argv[1]
    bampath2 = sys.argv[2]

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

    print("Perfect1: ", len(perfect1))
    print("Perfect2: ", len(perfect2))
    print("P1 & P2: ", len(perfect1 & perfect2))
    print("P1 - P2: ", len(perfect1 - perfect2))
    print("(P1-P2)/P1: ", round(len(perfect1 - perfect2)/len(perfect1)*100, 2))

if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "pre":
        precision()
    elif mode == "rec":
        recall()
    elif mode == "exunc":
        extract_uncovering()
