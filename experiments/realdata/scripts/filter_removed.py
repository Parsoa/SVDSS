import sys
from pysam import VariantFile


def parse_log(fpath):
    removed = set()
    for line in open(fpath):
        if not line.startswith("The site"):
            continue
        vinfo = line.split(' ')[2]
        removed.add(vinfo)
    return removed


def main():
    vcf_path = sys.argv[1]
    log1_path = sys.argv[2]
    log2_path = sys.argv[3]

    removed1 = parse_log(log1_path)
    removed2 = parse_log(log2_path)

    removed = removed1 | removed2

    vcf = VariantFile(vcf_path)
    print('\n'.join(str(vcf.header).split('\n')), end='')
    for record in vcf.fetch():
        chrom = record.contig
        pos = record.pos
        if f"{chrom}:{pos}" in removed:
            continue
        print(record, end='')

    print(len(removed1), len(removed2), len(
        removed1 & removed2), file=sys.stderr)


if __name__ == "__main__":
    main()
