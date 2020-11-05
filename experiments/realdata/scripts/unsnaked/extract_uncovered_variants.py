import sys

import re
from itertools import chain

from pysam import VariantFile

def main_output_vcf():
    bed1_path = sys.argv[1]
    bed2_path = sys.argv[2]
    vcf_path = sys.argv[3]

    idxs = {}
    for line in open(bed1_path):
        line = line.strip('\n').split('\t')
        idx, vlen, cov = int(line[3].split('_')[0]), len(line[5]), int(line[-1])
        if abs(vlen) > 1:
            idxs[idx] = False
            if cov > 0:
                idxs[idx] = idxs[idx] or True

    # for line in open(bed2_path):
    #     line = line.strip('\n').split('\t')
    #     idx, vlen, cov = line[3], len(line[5]), int(line[-1])
    #     if abs(vlen) > 1:
    #         if idx not in idxs:
    #             idxs[idx] = False
    #         if cov > 0:
    #             idxs[idx] = idxs[idx] or True

    vcf = VariantFile(vcf_path)
    i = 1
    print(vcf.header, end='')
    for record in vcf.fetch():
        if i in idxs and not idxs[i]:
            print(record, end='')
        i+=1

def main():
    bed_path = sys.argv[1]

    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        idx, vseq, cov = line[2], line[5], int(line[-1])
        vlen = len(vseq)
        if vlen > 1 and vlen < 250 and cov == 0:
            # check str
            pass

if __name__ == "__main__":
    main_output_vcf()
