import sys

from pysam import VariantFile

def main():
    vcf_path = sys.argv[1]
    sample1 = sys.argv[2]
    sample2 = sys.argv[3]
    out_prefix = sys.argv[4]

    out1_path = out_prefix + "_1.bed"
    out2_path = out_prefix + "_2.bed"

    out1 = open(out1_path, 'w')
    out2 = open(out2_path, 'w')

    vcf = VariantFile(vcf_path)
    idx = 1
    last_chrom = ''
    offset1 = 0
    offset2 = 0
    for record in vcf.fetch():
        chrom = record.contig
        if chrom != last_chrom:
            last_chrom = chrom
            offset1 = 0
            offset2 = 0

        pos = record.pos
        ref = record.ref
        alts = record.alts

        gt1, gt2 = record.samples[sample1]["GT"]
        gt1_2, gt2_2 = record.samples[sample2]["GT"]

        print1_flag = gt1 != gt1_2 and gt1 != gt2_2
        print2_flag = gt2 != gt1_2 and gt2 != gt2_2

        alt1 = ref
        if gt1 != 0:
            alt1 = alts[gt1-1]
        alt2 = ref
        if gt2 != 0:
            alt2 = alts[gt2-1]

        pos1 = pos-1 + offset1
        epos1 = pos1 + len(alt1)
        offset1 += len(alt1)-len(ref)
        if print1_flag:
            out1.write(f"{chrom}\t{pos1}\t{epos1}\t{idx}\t{ref}\t{alt1}\t{len(alt1)-len(ref)}\t{gt1}\n")

        pos2 = pos-1 + offset2
        epos2 = pos2 + len(alt2)
        offset2 += len(alt2)-len(ref)
        if print2_flag:
            out2.write(f"{chrom}\t{pos2}\t{epos2}\t{idx}\t{ref}\t{alt2}\t{len(alt2)-len(ref)}\t{gt2}\n")

        idx+=1

    out1.close()
    out2.close()
    
if __name__ == "__main__":
    main()
