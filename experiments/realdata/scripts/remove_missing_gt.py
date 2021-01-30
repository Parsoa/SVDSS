import sys
from pysam import VariantFile


def main():
    vcfpath = sys.argv[1]
    vcf = VariantFile(vcfpath)

    print('\n'.join(str(vcf.header).split('\n'))[:-1])
    hgnones = {0: 0, 1: 0, 2: 0}
    nanones = {0: 0, 1: 0, 2: 0}
    for record in vcf.fetch():
        hg1, hg2 = record.samples["HG00733"]["GT"]
        na1, na2 = record.samples["NA19240"]["GT"]
        hg = (hg1 is None) + (hg2 is None)
        na = (na1 is None) + (na2 is None)
        hgnones[hg] += 1
        nanones[na] += 1
        if hg1 is None or hg2 is None or na1 is None or na2 is None:
            continue

        print(record, end='')

    print(hgnones, file=sys.stderr)
    print(nanones, file=sys.stderr)


if __name__ == "__main__":
    main()
