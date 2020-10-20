import sys

import pysam

"""
"bedtools bamtobed" cannot set the quality to #errors+#deletions+#clips
"""

def get_errors(al, read_len):
    # bbmap:
    good_bases = al.get_cigar_stats()[0][7]
    # minimap2:
    # good_bases = al.get_cigar_stats()[0][0] - al.get_tag("NM")
    bad_bases = read_len - good_bases
    deletions = al.get_cigar_stats()[0][2]
    return bad_bases + deletions

def main():
    sampath = sys.argv[1]

    read_lens = {} # bbmap doesn't report read sequence on secondary alignment
    sam = pysam.AlignmentFile(sampath, 'r')
    skipped = 0
    for al in sam.fetch():
        ridx, ref, al_start, al_end = al.query_name, al.reference_name, al.reference_start, al.reference_end
        strand = '-' if al.is_reverse else '+'

        if ref is None:
            continue

        read_len = al.query_length
        if al.is_secondary:
            read_len = read_lens[al.query_name]
        else:
            read_lens[al.query_name] = read_len
        errors = get_errors(al, read_len)

        print(ref, al_start, al_end, ridx, errors, strand, sep='\t')
    print(f"Skipped {skipped} alignments (no primary)", file=sys.stderr)

if __name__ == "__main__":
    main()
