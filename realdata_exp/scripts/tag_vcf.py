import sys

from pysam import VariantFile

def classify_var(ref, alt):
    if ref == alt:
        return "HOMO"

    t = "SNP"
    if len(ref) != 1 or len(alt) != 1:
        if len(ref) < len(alt):
            t = "INS"
        elif len(ref) > len(alt):
            t = "DEL"
        else:
            t = "INDEL"
    return t

def tag():
    vcf_path = sys.argv[1]
    sample = sys.argv[2]

    vcf = VariantFile(vcf_path)

    vcf.header.add_line(f"##INFO=<ID=HAP1_{sample},Number=1,Type=String,Description=\"Type of variant on haplotype 1 ({sample})\">")
    vcf.header.add_line(f"##INFO=<ID=HAP2_{sample},Number=1,Type=String,Description=\"Type of variant on haplotype 2 ({sample})\">")

    print('\n'.join(str(vcf.header).split('\n'))[:-1])
    
    for record in vcf.fetch():
        ref = record.ref
        alts = record.alts
        gt1, gt2 = record.samples[sample]["GT"]

        alt1 = ref
        if gt1 != 0:
            alt1 = alts[gt1-1]
        alt2 = ref
        if gt2 != 0:
            alt2 = alts[gt2-1]

        htype1 = classify_var(ref, alt1)
        htype2 = classify_var(ref, alt2)

        record.info.__setitem__(f"HAP1_{sample}", htype1)
        record.info.__setitem__(f"HAP2_{sample}", htype2)

        print(record, end='')

def lens():
    vcf_path = sys.argv[1]
    sample = sys.argv[2]

    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        ref = record.ref
        alts = record.alts
        gt1, gt2 = record.samples[sample]["GT"]

        alt1 = ref
        if gt1 != 0:
            alt1 = alts[gt1-1]
        alt2 = ref
        if gt2 != 0:
            alt2 = alts[gt2-1]

        print(len(ref)-len(alt1), len(ref)-len(alt2), sep='\t')

def vcf2bed():
    vcf_path = sys.argv[1]
    sample = sys.argv[2]

    out1_path = vcf_path + "." + sample + ".hap1.bed"
    out2_path = vcf_path + "." + sample + ".hap2.bed"

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
            offset = 0

        pos = record.pos
        ref = record.ref
        alts = record.alts
        gt1, gt2 = record.samples[sample]["GT"]

        alt1 = ref
        if gt1 != 0:
            alt1 = alts[gt1-1]
        alt2 = ref
        if gt2 != 0:
            alt2 = alts[gt2-1]

        pos1 = pos-1 + offset1
        epos1 = pos1 + len(alt1)
        offset1 += len(alt1)-len(ref)
        out1.write(f"{chrom}_1\t{pos1}\t{epos1}\t{idx}\t{ref}\t{alt1}\t{len(alt1)-len(ref)}\t{gt1}\n")

        pos2 = pos-1 + offset1
        epos2 = pos2 + len(alt2)
        offset2 += len(alt2)-len(ref)
        out2.write(f"{chrom}_2\t{pos2}\t{epos2}\t{idx}\t{ref}\t{alt2}\t{len(alt2)-len(ref)}\t{gt2}\n")

        idx+=1

    out1.close()
    out2.close()

def vcf2bed_unique():
    vcf_path = sys.argv[1]
    sample1 = sys.argv[2]
    sample2 = sys.argv[3]

    out1_path = vcf_path + "." + sample1 + ".unique_wrt_" + sample2 + ".hap1.bed"
    out2_path = vcf_path + "." + sample1 + ".unique_wrt_" + sample2 + ".hap2.bed"

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
    
def extract_unique_vcf():
    vcf_path = sys.argv[1]
    sample1 = sys.argv[2]
    sample2 = sys.argv[3]

    vcf = VariantFile(vcf_path)
    print('\n'.join(str(vcf.header).split('\n'))[:-1])
    for record in vcf.fetch():
        gt1_1, gt2_1 = record.samples[sample1]["GT"]
        gt1_2, gt2_2 = record.samples[sample2]["GT"]
        if (gt1_1 != gt1_2 and gt1_1 != gt2_2) or (gt2_1 != gt1_2 and gt2_1 != gt2_2):
            print(record, end='')
    
if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "unq":
        extract_unique_vcf()
    elif mode == "v2b":
        vcf2bed()
    elif mode == "v2bu":
        vcf2bed_unique()
