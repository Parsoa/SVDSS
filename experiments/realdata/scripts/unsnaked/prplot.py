import sys
import copy

import matplotlib.pyplot as plt

def get_abundance(ridx):
    return int(ridx.split('#')[-1])

def extract_truth(bed_path):
    truth = set()
    for line in open(bed_path, 'r'):
        truth.add(int(line.split('\t')[3]))
    return truth

def extract_variants(bed_path):
    variants = {}
    for line in open(bed_path, 'r'):
        line = line.strip('\n').split('\t')
        vidx, ridx = int(line[3]), line[11]
        if vidx not in variants:
            variants[vidx] = set()
        variants[vidx].add(ridx)
    return variants

def combine_variants(variants1, variants2):
    variants = copy.deepcopy(variants1)
    for vidx,ridxs in variants2.items():
        if vidx in variants:
            variants[vidx] |= ridxs
        else:
            variants[vidx] = ridxs
    return variants

def compute_recall(truth, variants, minab=0):
    found = set()
    for vidx, ridxs in variants.items():
        maxab = max([get_abundance(ridx) for ridx in ridxs])
        if maxab < minab:
            continue
        found.add(vidx)

    total = len(truth)
    R = round(len(found)/total*100, 2) if total != 0 else 0
    print(minab, len(found), total, R, "(R)")
    return R

def extract_alignments(bed_path):
    alignments = {}
    for line in open(bed_path, 'r'):
        line = line.strip('\n').split('\t')
        ridx, cov = line[3], int(line[-1])
        if ridx not in alignments:
            alignments[ridx] = False
        alignments[ridx] = alignments[ridx] or cov != 0
    return alignments

def combine_alignments(alignments1, alignments2):
    alignments = copy.deepcopy(alignments1)

    for ridx, cov_flag in alignments2.items():
        if ridx in alignments:
            alignments[ridx] = alignments[ridx] or cov_flag
        else:
            alignments[ridx] = cov_flag2
    return alignments

def compute_precision(alignments, minab = 0):
    total = 0
    covering = 0
    for ridx, cov_flag in alignments.items():
        if get_abundance(ridx) < minab:
            continue
        total += 1
        covering += cov_flag
    P = round(covering/total*100, 2) if total > 0 else 0
    print(minab, covering, total, P, "(P)")
    return P

def main():
    bed_truth1_path = sys.argv[1]
    bed_truth2_path = sys.argv[2]
    bed_rec1_path = sys.argv[3]
    bed_rec2_path = sys.argv[4]
    bed_pre1_path = sys.argv[5]
    bed_pre2_path = sys.argv[6]

    truth1 = extract_truth(bed_truth1_path)
    truth2 = extract_truth(bed_truth2_path)
    truth = truth1 | truth2

    variants1 = extract_variants(bed_rec1_path)
    variants2 = extract_variants(bed_rec2_path)
    variants = combine_variants(variants1, variants2)

    alignments1 = extract_alignments(bed_pre1_path)
    alignments2 = extract_alignments(bed_pre2_path)
    alignments = combine_alignments(alignments1, alignments2)

    for i in range(5,9):        
        P = compute_precision(alignments, i)
        # Recall by covered variants:
        # R = compute_recall(truth, variants, i)

        # Recalls by covered alleles:
        R1 = compute_recall(truth1, variants1, i)
        R2 = compute_recall(truth1, variants2, i)
        print("")

def plot():
    Xs = range(5,18)
    # Recall by variant results (run main to get these values)
    Ps = [62.52, 63.47, 63.63, 63.33, 62.73, 61.96, 61.14, 60.29, 59.34, 58.2, 56.62, 54.44, 51.47] #, 47.61, 42.81, 37.18, 31.06, 24.81, 19.03, 14.21, 10.31, 7.53, 5.56, 4.25, 3.49, 3.06]
    Rs = [93.22, 88.49, 80.96, 71.12, 60.27, 50.03, 41.35, 34.44, 28.98, 24.42, 20.34, 16.56, 13.07] #, 9.98, 7.33, 5.2, 3.57, 2.38, 1.56, 1.02, 0.66, 0.44, 0.3, 0.21, 0.16, 0.14]

    fig, ax = plt.subplots(1, 1, tight_layout=True)
    ax.scatter(Ps, Rs, s=5)

    for i, txt in enumerate(Xs):
        ax.annotate(txt, (Ps[i], Rs[i]))
    ax.set_xlabel("Precision")
    ax.set_ylabel("Recall")
    ax.set_xlim(0,100)
    ax.set_ylim(0,100)
    plt.show()

if __name__ == "__main__":
    main()
