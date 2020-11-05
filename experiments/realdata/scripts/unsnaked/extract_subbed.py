import sys

def main():
    bedpath = sys.argv[1]
    fpath = sys.argv[2]

    idxs = set()
    for line in open(fpath):
        idxs.add(line.strip('\n'))

    for line in open(bedpath):
        idx = line.strip('\n').split('\t')[3]
        if idx in idxs:
            print(line, end='')

if __name__ == "__main__":
    main()
