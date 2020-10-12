import sys

from Bio import SeqIO

def main():
    fa_path = sys.argv[1]
    min_l = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    L = 0
    N = 0
    for record in SeqIO.parse(fa_path, "fasta"):
        seq = str(record.seq).upper()
        l = len(seq)
        if l >= min_l:
            n = seq.count('N')
            L += l
            N += n
    print(N, L, N/L)

if __name__ == "__main__":
    main()
