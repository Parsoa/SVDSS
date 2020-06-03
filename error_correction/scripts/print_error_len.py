import sys, os
import gzip

from Bio import SeqIO

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
    last_idx = ""
    last_end = 0
    X = []
    for record in parse_file(fpath):
        idx = record.id.strip('\n').split('/')[:-1]
        s,e = record.id.strip('\n').split('_')[-1].split(':')
        s,e = int(s),int(e)

        if idx != last_idx:
            last_idx = idx
        else:
            X.append(s-last_end)
            print(s-last_end)
        last_end = e

if __name__ == "__main__":
    main()
