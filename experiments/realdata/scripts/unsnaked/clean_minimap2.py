import sys

import pysam

def main():
    sampath = sys.argv[1]

    sam = pysam.AlignmentFile(sampath, 'r')
    print(sam.header, end='')

    alignments = set()
    for al in sam.fetch():
        if al.query_name not in alignments:
            alignments.add(al.query_name)
            print(al.tostring(sam))
        else:
            if al.is_secondary or al.is_supplementary or al.is_unmapped:
                print(al.tostring(sam))
            else:
                al.flag = 272 if al.is_reverse else 256
                print(al.tostring(sam))
    sam.close()

if __name__ == "__main__":
    main()
