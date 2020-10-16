import sys

from Bio import SeqIO
import pysam

def main():
    fa1_path = sys.argv[1] # contigs hap1
    fa2_path = sys.argv[2] # contigs hap2
    sampath = sys.argv[3]  # alignments
    nal = 2

    Idxs1 = set()
    for record in SeqIO.parse(fa1_path, "fasta"):
        Idxs1.add(record.id)

    Idxs2 = set()
    for record in SeqIO.parse(fa2_path, "fasta"):
        idx = record.id
        if idx in Idxs1:
            idx += "_2"
        Idxs2.add(idx)
    
    aligns = {}
    sam = pysam.AlignmentFile(sampath, "r")
    for al in sam.fetch():
        if al.query_name not in aligns:
            aligns[al.query_name] = []
        aligns[al.query_name].append(al.reference_name)

    total = 0
    diff_contigs = 0
    diff_contigs_same_cluster = 0
    diff_contigs_same_cluster_same_fragment = 0
    for qidx, ridxs in aligns.items():
        if len(ridxs) == nal:
            total += 1
            idx1, idx2 = ridxs
            if (idx1 in Idxs1 and idx2 in Idxs2) or (idx1 in Idxs2 and idx2 in Idxs1):
                diff_contigs += 1
                cl1, cl2 = idx1.split('_')[0], idx2.split('_')[0]
                if cl1 == cl2:
                    diff_contigs_same_cluster += 1
                    f1, f2 = idx1.split('_')[1], idx2.split('_')[1]
                    if f1 == f2:
                        diff_contigs_same_cluster_same_fragment += 1

    print("Diff contigs:", diff_contigs, total, diff_contigs/total)
    print("Diff contigs, same cluster:", diff_contigs_same_cluster, total, diff_contigs_same_cluster/total)
    print("Diff contigs, same cluster, same fragment", diff_contigs_same_cluster_same_fragment, total, diff_contigs_same_cluster_same_fragment/total)

if __name__ == "__main__":
    main()
