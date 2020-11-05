import sys

import pysam

def main():
    sampath = sys.argv[1]
    nal = 5 # int(sys.argv[2])

    aligns = {}
    quals = {}
    seqs = {}
    sam = pysam.AlignmentFile(sampath, "r")
    for al in sam.fetch():
        aligns[al.query_name] = aligns[al.query_name]+1 if al.query_name in aligns else 1
        if not al.is_secondary:
            good_bases = al.get_cigar_stats()[0][7]
            bad_bases = al.query_length - good_bases
            deletions = al.get_cigar_stats()[0][2]
            quals[al.query_name] = bad_bases + deletions
            seqs[al.query_name] = al.query_sequence

    for idx,count in aligns.items():
        if count == nal:
            print(f">{idx}_{quals[idx]}\n{seqs[idx]}")

if __name__ == "__main__":
    main()
