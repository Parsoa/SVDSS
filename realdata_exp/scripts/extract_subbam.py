import sys

import pysam

def main():
    fpath = sys.argv[1]
    bampath = sys.argv[2]

    idxs = set()
    for line in open(fpath):
        idxs.add(line.strip('\n'))

    bam = pysam.AlignmentFile(bampath, "r")
    print(bam.header, end='')
    for record in bam.fetch():
        if record.query_name in idxs:
            print(record.tostring(bam))

if __name__ == "__main__":
    main()
