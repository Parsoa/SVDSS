import sys

from Bio.Seq import Seq
import pysam

def main():
    bampath = sys.argv[1]

    nerr = sys.argv[2]
    exact = True
    if nerr.startswith('>'):
        exact = False
        nerr = nerr[1:]
    nerr = int(nerr)

    out_fmt = sys.argv[3] if len(sys.argv) == 4 else "sam"

    good_als = 0
    bam = pysam.AlignmentFile(bampath, "r")
    if out_fmt == "sam":
        print(bam.header, end='')
    for al in bam.fetch():
        good_bases = al.get_cigar_stats()[0][7]
        bad_bases = al.query_length - good_bases
        deletions = al.get_cigar_stats()[0][2]
        if (exact and bad_bases + deletions == nerr) or (not exact and bad_bases + deletions > nerr):
            good_als += 1
            if out_fmt == "sam":
                print(al.tostring(bam))
            elif out_fmt == "fa":
                seq = Seq(al.query_sequence)
                if al.is_reverse:
                    seq = seq.reverse_complement()
                print(f">{al.query_name}\n{seq}")

    print("Total good alignments:", good_als, file=sys.stderr)

if __name__ == "__main__":
    main()
