import sys

import pysam

def main():
    sampath = sys.argv[1]
    fpath = sys.argv[2]

    idxs = set()
    for line in open(fpath):
        idxs.add(line.strip('\n'))

    sam = pysam.AlignmentFile(sampath, 'r')
    print(sam.header, end='')
    for al in sam.fetch():
        if al.query_name in idxs:
            print(al.tostring(sam))

if __name__ == "__main__":
    main()
