import sys

def main():
    fpath = sys.argv[1]
    l_thresh = int(sys.argv[2])

    found_small = 0
    found_big = 0
    tot_small = 0
    tot_big = 0
    for line in open(fpath):
        idx, l, ov = line.strip('\n').split(' ')
        l, ov = abs(int(l)), int(ov)

        if l < l_thresh:
            tot_small += 1
            if ov > 0:
                found_small += 1
        else:
            tot_big += 1
            if ov > 0:
                found_big += 1

    print(f"Split at", l_thresh)
    print("Type", "Found", "Total", "%", sep='\t')
    print("<", found_small, tot_small, found_small/tot_small, sep='\t')
    print(">=", found_big, tot_big, found_big/tot_big, sep='\t')
    print("All", found_small + found_big, tot_small + tot_big, (found_small + found_big)/(tot_small + tot_big), sep='\t')

if __name__ == "__main__":
    main()
