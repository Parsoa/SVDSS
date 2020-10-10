import sys

import pysam

def main():
    bampath = sys.argv[1]
    fpath = sys.argv[2]

    idxs = set()
    for line in open(fpath):
        idxs.add(line.strip('\n'))

    bam = pysam.AlignmentFile(bampath, "rb")
    print(bam.header, end='')
    for record in bam.fetch():
        if record.query_name in idxs:
            print(record.tostring(bam))

if __name__ == "__main__":
    main()
