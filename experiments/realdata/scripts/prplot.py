import sys

from pysam import VariantFile

"""
mkdir prplot ; cd prplot

bedtools intersect -wao -a $data/pingpong_ext/OUT/HGspecificvars_1.bed -b $data/pingpong_ext/OUT/pp/over/sspecific.short.hg-haps.sam.bed > rec1.short.bed
bedtools intersect -wao -a $data/pingpong_ext/OUT/HGspecificvars_2.bed -b $data/pingpong_ext/OUT/pp/over/sspecific.short.hg-haps.sam.bed > rec2.short.bed
bedtools intersect -wao -a $data/pingpong_ext/OUT/HGspecificvars_1.bed -b $data/pingpong_ext/OUT/pp/over/sspecific.long.hg-haps.sam.bed > rec1.long.bed
bedtools intersect -wao -a $data/pingpong_ext/OUT/HGspecificvars_2.bed -b $data/pingpong_ext/OUT/pp/over/sspecific.long.hg-haps.sam.bed > rec2.long.bed
cat rec1.short.bed rec1.long.bed > rec1.bed
cat rec2.short.bed rec2.long.bed > rec2.bed

bedtools intersect -wao -b $data/pingpong_ext/OUT/HGspecificvars_1.bed -a $data/pingpong_ext/OUT/pp/over/sspecific.short.hg-haps.sam.bed > pre1.short.bed
bedtools intersect -wao -b $data/pingpong_ext/OUT/HGspecificvars_2.bed -a $data/pingpong_ext/OUT/pp/over/sspecific.short.hg-haps.sam.bed > pre2.short.bed
bedtools intersect -wao -b $data/pingpong_ext/OUT/HGspecificvars_1.bed -a $data/pingpong_ext/OUT/pp/over/sspecific.long.hg-haps.sam.bed > pre1.long.bed
bedtools intersect -wao -b $data/pingpong_ext/OUT/HGspecificvars_2.bed -a $data/pingpong_ext/OUT/pp/over/sspecific.long.hg-haps.sam.bed > pre2.long.bed
cat pre1.short.bed pre1.long.bed > pre1.bed
cat pre2.short.bed pre2.long.bed > pre2.bed
"""

def parse_recbed(bed_path):
    alleles = {}
    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        ref, vidx, ridx, covering = line[0], line[3], line[8], int(line[-1])
        aidx = vidx + ('2' if ref.endswith("_2") else '1')
        if aidx not in alleles:
            alleles[aidx] = set()
        if covering > 0:
            alleles[aidx].add(ridx)
    return alleles

def recall(alleles, cutoff):
    total = 0
    hits = 0
    for ridxs in alleles.values():
        total += 1
        if any([int(ridx.split('#')[1]) >= cutoff for ridx in ridxs]):
            hits += 1
    return hits, total

# precision setup
def count_sample(fq_path):
    # assuming fastq
    n = 4
    reads_by_coverage = {}
    for line in open(fq_path):
        if n % 4 == 0:
            line = line[1:-1] # removes @ and \n
            cov = int(line.split('#')[1])
            reads_by_coverage[cov] = reads_by_coverage[cov] + 1 if cov in reads_by_coverage else 1
            n = 0
        n += 1
    return reads_by_coverage

def parse_vcf(vcf_path, sample1="HG00733", sample2="NA19240"):
    variants = {}
    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        idx = record.id
        gt1, gt2 = record.samples[sample1]["GT"]
        gt1_2, gt2_2 = record.samples[sample2]["GT"]
        variants[idx] = [(gt1, gt2), (gt1_2, gt2_2)]
    return variants

def parse_prebed(bed_path, firsthap):
    alignments = {}
    for line in open(bed_path):
        line = line.strip('\n').split('\t')
        ref, start, end, ridx, vidx, covering = line[0], line[1], line[2], line[3], line[9], int(line[-1])
        if ref.endswith("_2") and firsthap:
            continue
        posidx = ridx + ":" + start + "-" + end
        if posidx not in alignments:
            alignments[posidx] = set()
        if covering > 0:
            if vidx == ".":
                print(line)
                exit(1)
            alignments[posidx].add(vidx)
    return alignments

def check_haplo_uniqueness(variants, firsthap=True):
    HGhaplo = [v[0][not firsthap] for v in variants]
    NAhaplo1 = [v[1][0] for v in variants]
    NAhaplo2 = [v[1][1] for v in variants]
    return HGhaplo != NAhaplo1 and HGhaplo != NAhaplo2

def precision(alignments1, alignments2, variants, cutoff):
    covering_results = {}
    for posidx, vidxs in alignments1.items():
        ridx = posidx.split(':')[0]
        if int(ridx.split('#')[1]) < cutoff:
            continue
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, True)
            covering_results[ridx] = covering_results[ridx] or unique

    for posidx, vidxs in alignments2.items():
        ridx = posidx.split(':')[0]
        if int(ridx.split('#')[1]) < cutoff:
            continue
        if ridx not in covering_results:
            covering_results[ridx] = False
        if len(vidxs) > 0:
            covered_variants = [variants[vidx] for vidx in vidxs]
            unique = check_haplo_uniqueness(covered_variants, False)
            covering_results[ridx] = covering_results[ridx] or unique

    covering = 0
    for ridx, flag in covering_results.items():
        if flag:
            covering += 1
    return covering

def main():
    fq_path = sys.argv[1]
    vcf_path = sys.argv[2]
    rec1_waobed = sys.argv[3] # bedtools intersect -wao {variants} {alignments}
    rec2_waobed = sys.argv[4]
    pre1_waobed = sys.argv[5] # bedtools intersect -wao {alignments} {variants}
    pre2_waobed = sys.argv[6]

    # # RECALL
    # alleles1 = parse_recbed(rec1_waobed)
    # print("Parsed rec1", file=sys.stderr)
    # alleles2 = parse_recbed(rec2_waobed)
    # print("Parsed rec2", file=sys.stderr)
    # alleles = {**alleles1, **alleles2}
    # print("Combined", file=sys.stderr)
    # for c in range(5, 11):
    #     print(f"Computing recall with c={c}", file=sys.stderr)
    #     hits, total = recall(alleles, c)
    #     R = round(hits/total*100 if total > 0 else 0, 3)
    #     print("R", c, hits, total, R, sep=',')

    # PRECISION
    print("Parsing sample", file=sys.stderr)
    reads_by_coverage = count_sample(fq_path)
    print("Parsing VCF", file=sys.stderr)
    variants = parse_vcf(vcf_path)
    print("Parsing pre1", file=sys.stderr)
    alignments1 = parse_prebed(pre1_waobed, True)
    print("Parsing pre2", file=sys.stderr)
    alignments2 = parse_prebed(pre2_waobed, False)

    for c in range(5, 11):
        print(f"Computing recall with c={c}", file=sys.stderr)
        covering = precision(alignments1, alignments2, variants, c)

        total = 0 # TODO improve this
        for cutoff,n in reads_by_coverage.items():
            if cutoff >= c:
                total += n

        P = round(covering/total*100 if total != 0 else 0, 3)
        print("P", c, covering, total, P, sep=',')

if __name__ == "__main__":
    main()
