import sys

from Bio import SeqIO

def main():
    fqpath = sys.argv[1]
    fpath = sys.argv[2]

    idxs = set()
    for line in open(fpath):
        idxs.add(line.strip('\n'))

    for record in SeqIO.parse(fqpath, "fastq"):
        if record.id in idxs:
            SeqIO.write(record, sys.stdout, "fastq")

if __name__ == "__main__":
    main()
