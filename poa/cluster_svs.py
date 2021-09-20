import sys
from pysam import VariantFile
from intervaltree import IntervalTree

def cluster_by_len(variations):
    l_clusters = {}
    for v in variations:
        new_flag = True
        for l,vs in l_clusters.items():
            if abs(v.info["SVLEN"] - l) <= 10:
                vs.append(v)
                new_flag = False
                break
        if new_flag:
            l_clusters[v.info["SVLEN"]] = [v]
    merged_vs = []
    for _,vs in l_clusters.items():
        new_w = 0
        for v in vs:
            new_w += v.info["WEIGHT"]
        vs[0].info["WEIGHT"] = new_w
        merged_vs.append(vs[0])
    return merged_vs

def main():
    vcf_path = sys.argv[1]
    radius = 250

    tree = IntervalTree()
    vcf = VariantFile(vcf_path)
    for rec in vcf.fetch():
        s = rec.pos
        e = rec.stop + 1
        overlaps = tree.overlap(s, e)
        if len(overlaps) == 0:
            tree.addi(s - radius, e + radius, [rec])
        else:
            minb = s - radius
            maxe = e + radius
            data = [rec]
            for overlap in overlaps:
                minb = min(minb, overlap.begin)
                maxe = max(maxe, overlap.end)
                data.extend(overlap.data)
                tree.remove(overlap)
            tree.addi(minb, maxe, data)

    print(vcf.header, end='')
    for interval in tree:
        if len(interval.data) == 1:
            print(interval.data[0], end='')
        else:
            insertions = []
            deletions = []
            for rec in interval.data:
                if rec.info["SVLEN"] < 0:
                    deletions.append(rec)
                else:
                    insertions.append(rec)

            # DELETIONS
            if len(deletions) == 1:
                print(deletions[0], end='')
            else:
                # same size -> same variation
                deletions = cluster_by_len(deletions)

                # same ratio -> same variation
                for d in deletions:
                    print(d, end='')
            
            if len(deletions) > 2:
                print(deletions[0].pos, end=' ', file=sys.stderr)
                for d in deletions:
                    print(d.info["SVLEN"], d.info["WEIGHT"]/d.info["COV"], end=' ', file=sys.stderr)
                print('', file=sys.stderr)
            
            # INSERTIONS
            if len(insertions) == 1:
                print(insertions[0], end='')
            else:
                # same size -> same variation
                insertions = cluster_by_len(insertions)
                for i in insertions:
                    print(i, end='')

            if len(insertions) > 2:
                print(insertions[0].pos, end=' ', file=sys.stderr)
                for i in insertions:
                    print(i.info["SVLEN"], i.info["WEIGHT"]/i.info["COV"], end=' ', file=sys.stderr)
                print('', file=sys.stderr)
            # same w -> combine
            # diff w -> merge

if __name__ == "__main__":
    main()