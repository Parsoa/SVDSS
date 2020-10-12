import sys

def main():
    fpath = sys.argv[1]
    l_thresh = int(sys.argv[2])

    found = {"SNP" : 0, "<" : 0, ">=" : 0}
    tot = {"SNP" : 0, "<" : 0, ">=" : 0}
    for line in open(fpath):
        idx, l, ov = line.strip('\n').split(' ')
        l, ov = abs(int(l)), int(ov)

        idx = "SNP"
        if l == 0:
            idx = "SNP"
        elif l < l_thresh:
            idx = "<"
        else:
            idx = ">="

        tot[idx] += 1
        if ov > 0:
            found[idx] += 1

    print(f"Split at", l_thresh)
    print("")
    print("Type", "Found", "Total", "%", sep='\t')
    print("---")
    for idx in ["SNP", "<", ">="]:
        print(idx, found[idx], tot[idx], round(found[idx]/tot[idx]*100, 2), sep='\t')
    print("---")
    print("All", sum(found.values()), sum(tot.values()), round(sum(found.values())/sum(tot.values())*100, 2), sep='\t')

if __name__ == "__main__":
    main()
