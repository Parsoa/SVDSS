import sys

import pysam

"""
"bedtools bamtobed" cannot set the quality to #errors+#clips (it's more convenient for analysis)
"""

def main():
    bam_path = sys.argv[1]

    bam = pysam.AlignmentFile(bam_path, "rb")
    for al in bam.fetch():
        ridx, ref, al_start, al_end = al.query_name, al.reference_name, al.reference_start, al.reference_end
        strand = '-' if al.is_reverse else '+'

        if ref is None:
            continue

        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases

        print(ref, al_start, al_end, ridx, bad_bases, strand, sep='\t')

if __name__ == "__main__":
    main()
