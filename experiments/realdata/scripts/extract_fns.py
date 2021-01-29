import sys

from Bio import SeqIO
import pysam
from pysam import VariantFile

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

###########################
### ALLELE-BASED RECALL ###
###########################


def get_vtype(idx, splitL):
    idx = idx.split('-')
    if idx[2] == "SNV":
        return "SNP"
    else:
        if int(idx[-1]) < splitL:
            return "INDEL"
        else:
            return "SV"


def extract_alleles(bed_path, splitL):
    alleles = {}
    for line in open(bed_path, 'r'):
        line = line.strip('\n').split('\t')
        chrom, idx, gt, covering = line[0], line[3], int(
            line[-2]), int(line[-1])
        vtype = get_vtype(idx, splitL)
        hap = '2' if chrom.endswith("_2") else '1'
        idx += hap
        gt = "Ref" if gt == 0 else "Alt"
        if idx in alleles:
            alleles[idx] = (covering > 0 or alleles[idx][0], vtype, hap, gt)
        else:
            alleles[idx] = (covering > 0, vtype, hap, gt)
    return alleles


def main():
    # TODO this can be done better (more pythonic)
    # bedintersect -c between alignments and unique variants on first haplotype
    bed1_path = sys.argv[1]
    # bedintersect -c between alignments and unique variants on second haplotype
    bed2_path = sys.argv[2]
    vcf_path = sys.argv[3]
    vtype = sys.argv[4]
    splitL = 50

    alleles1 = extract_alleles(bed1_path, splitL)
    alleles2 = extract_alleles(bed2_path, splitL)

    alleles = {**alleles1, **alleles2}

    fns = set()
    for aidx, (covering, vt, hap, gt) in alleles.items():
        if vt == vtype:
            if covering:
                fns.discard(aidx[:-1])
            else:
                fns.add(aidx[:-1])

    vcf = VariantFile(vcf_path)
    print(vcf.header, end='')
    for record in vcf.fetch():
        idx = record.id
        if idx in fns:
            print(record, end='')

if __name__ == "__main__":
    main()
