import sys
import copy

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def extract_from_bed(bed_path):
    variants = {}
    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        ref = line[0]
        if not (ref == "chr22" or ref == "chr22_2"):
            continue
        vidx, l, ridx = line[3], int(line[6]), line[-3]
        if vidx not in variants:
            variants[vidx] = set()
        variants[vidx].add(ridx)
    return variants

def combine_variants(variants1, variants2):
    variants = copy.deepcopy(variants1)
    for vidx,ridxs in variants2.items():
        if vidx in variants:
            variants[vidx] |= ridxs
        else:
            variants[vidx] = ridxs
    return variants

def get_data(variants):
    data = []
    for vidx,ridxs in variants.items():
        maxab = max([int(ridx.split('#')[-1]) for ridx in ridxs])
        data.append(maxab)
    return data

def main():
    bed1_path = sys.argv[1]
    bed2_path = sys.argv[2]

    variants1 = extract_from_bed(bed1_path)
    variants2 = extract_from_bed(bed2_path)
    variants = combine_variants(variants1, variants2)

    # for vidx,ridxs in variants.items():
    #     if len(ridxs) == 1 and int(list(ridxs)[0].split('#')[-1]) == 5:
    #         print(vidx, list(ridxs)[0])
    # return

    data1 = get_data(variants1)
    data2 = get_data(variants2)
    data12 = get_data(variants)

    ddata1 = {}
    for l in data1:
        ddata1[l] = ddata1[l]+1 if l in ddata1 else 1
    ddata2 = {}
    for l in data2:
        ddata2[l] = ddata2[l]+1 if l in ddata2 else 1
    ddata12 = {}
    for l in data12:
        ddata12[l] = ddata12[l]+1 if l in ddata12 else 1

    print('\t'.join(str(x) for x in range(5,31)))
    print('\t'.join([str(ddata1[l]) if l in ddata1 else "0" for l in range(5,31)]))
    print('\t'.join([str(ddata2[l]) if l in ddata2 else "0" for l in range(5,31)]))
    print('\t'.join([str(ddata12[l]) if l in ddata12 else "0" for l in range(5,31)]))

    data = {"Hap" : ['1']*len(data1) + ['2']*len(data2) + ["12"]*len(data12),
            "Abs" : data1 + data2 + data12}

    sns.histplot(data, x="Abs", hue="Hap",
                 bins = range(5,31),
                 element="step",
                 binrange = (5,30))

    plt.show()

if __name__ == "__main__":
    main()
