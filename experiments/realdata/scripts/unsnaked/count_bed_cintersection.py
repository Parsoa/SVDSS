import sys

def main():
    data = {}
    for line in open(sys.argv[1]):
        line = line.strip('\n').split('\t')
        idx, found = line[3], line[-1]
        found = int(found)
        if idx not in data:
            data[idx] = 0
        data[idx] = max(data[idx], found)

    total = len(data)
    over = 0
    for v in data.values():
        if v>0:
            over+=1
    print(over, total, over/total)

if __name__ == "__main__":
    main()
