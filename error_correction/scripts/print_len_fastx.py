import sys, os
import gzip

from Bio import SeqIO

FA  = 0
FAZ = 1
FQ  = 2
FQZ = 3

def parse_file(fpath):
    prev = []
    if fpath.strip('\n').split('.')[-1] == "gz":
        f = gzip.open(fpath, 'rt')
        prev = fpath.strip('\n').split('.')[-2]
    else:
        f = open(fpath, 'r')
        prev = fpath.strip('\n').split('.')[-1]

    if prev in ["fa", "fasta"]:
        return SeqIO.parse(f, "fasta")
    elif prev in ["fq", "fastq"]:
        return SeqIO.parse(f, "fastq")
    else:
        print("File extension not valid", file=sys.stderr)
        sys.exit(1)

def main():
    fpath = sys.argv[1]
    for record in parse_file(fpath):
        print(len(record))

if __name__ == "__main__":
    main()
