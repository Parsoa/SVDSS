import sys
import time
from intervaltree import IntervalTree
from Bio import SeqIO
import pysam


def getTime():
    return time.strftime('[ %b %d, %Y - %l:%M:%S%p ]')


def eprint(*args, **kwargs):
    print(getTime(), *args, file=sys.stderr, flush=True, **kwargs)


class SFS:
    def __init__(self, chrom, ridx, rs, re, qs, qe, seq):
        self.chrom = chrom  # chromosome name
        self.ridx = ridx  # read identifier
        self.rs = rs  # start on reference
        self.re = re  # end on reference
        self.qs = qs  # start on read (query)
        self.qe = qe  # end on read
        self.seq = seq  # sequence on read

    def __repr__(self):
        return f"{self.chrom}:{self.rs+1}-{self.re+1} {self.ridx} {self.seq}"

    def __hash__(self):
        return hash(f"{self.ridx}{self.rs}{self.re}")


class Cluster:
    def __init__(self, chrom):
        self.chrom = chrom
        self.s = float("inf")
        self.e = 0
        self.sfss = set()
        self.ridxs = set()
        self.cov = 0

    def add(self, sfs):
        if sfs.ridx not in self.ridxs:
            self.sfss.add(sfs)
            self.ridxs.add(sfs.ridx)
        if sfs.rs < self.s:
            self.s = sfs.rs
        if sfs.re > self.e:
            self.e = sfs.re

    def idx(self):
        return f"{self.chrom}:{self.s+1}-{self.e+1}"

    def get_seqs(self):
        return [sfs.seq for sfs in self.sfss]

    def set_coverage(self, cov):
        self.cov = cov

    def dump(self):
        print(self.chrom, self.s, self.e, self.cov, len(
            self.sfss), " ".join(self.get_seqs()), sep=" ")

    def __repr__(self):
        return f"{self.chrom}:{self.s+1}-{self.e+1} {len(self)}"

    def __len__(self):
        return len(self.sfss)

    def __iter__(self):
        return iter(self.sfss)


def get_unique_kmer(alpairs, qseq, rseq, k, fromend=True):
    kmers = {}
    Is = list(range(0, len(alpairs) - k + 1))
    for i in Is:
        kmer = alpairs[i:i+k]
        if any([q is None or r is None for (q, r) in kmer]):
            # we want clean kmers only - ie placed kmers, no insertions or deletions
            continue
        # qkmer_seq = qseq[kmer[0][0]:kmer[-1][0]+1]
        rkmer_seq = rseq[kmer[0][1]:kmer[-1][1]+1]
        # note: after read reconstruction these kmers should be the same

        kmers[rkmer_seq] = kmers[rkmer_seq] + 1 if rkmer_seq in kmers else 1

    if fromend:
        Is = Is[::-1]
    last_kmer = [(None, None)]
    last_rkmer_seq = ""
    for i in Is:
        kmer = alpairs[i:i+k]
        if any([q is None or r is None for (q, r) in kmer]):
            # we want clean kmers only - ie placed kmers, no insertions or deletions
            continue
        last_kmer = kmer
        rkmer_seq = rseq[kmer[0][1]:kmer[-1][1]+1]
        last_rkmer_seq = rkmer_seq
        if kmers[rkmer_seq] == 1:
            return kmer[0], rkmer_seq
    return last_kmer[0], last_rkmer_seq


def main():
    fa_path = sys.argv[1]
    sfs_path = sys.argv[2]
    bam_path = sys.argv[3]

    # some hardcoded parameters FIXME
    k = 7
    maxw = 100
    minsupp = 2

    # some stats on clips and clusters
    clips_1 = 0
    clips_2 = 0
    clips_3 = 0
    small_cl = 0
    small_extcl = 0
    extcl = 0

    # Reference
    eprint("Parsing reference..")
    chroms = {}
    for rec in SeqIO.parse(fa_path, "fasta"):
        chroms[rec.id] = str(rec.seq)

    # SFSs
    eprint("Parsing SFSs..")
    sfss = {}
    last_ridx = ""
    for line in open(sfs_path):
        ridx, p, l, c = line.strip("\n").split()
        if ridx != "*":
            sfss[ridx] = set()
            last_ridx = ridx
        p, l = int(p), int(l)
        sfss[last_ridx].add((p, p+l-1))

    # Get extended superstring from read alignments
    eprint("Parsing alignments..")
    extsupersfss = set()
    bam = pysam.AlignmentFile(bam_path, "rb")
    n_al = 0
    for aln in bam.fetch():
        if aln.is_supplementary or aln.is_unmapped or aln.is_secondary:
            continue
        qname = aln.query_name
        if qname not in sfss:
            continue
        chrom = aln.reference_name

        # pos = aln.reference_start
        qseq = aln.query_sequence
        alpairs = aln.get_aligned_pairs()

        for sfs in sfss[qname]:
            s, e = sfs[0], sfs[1]
            # 1 we get the first bases on each direction (downstream and upstream the supersting) that have been correctly placed
            refs, refe = -1, -1
            for q, r in alpairs[::-1]:
                if q == None or r == None:
                    continue
                if q <= s:
                    refs = r
                    break
            for q, r in alpairs:
                if q == None or r == None:
                    continue
                if q >= e:
                    refe = r
                    break

            # 2 we extract the local alignment of the region of interest
            local_alpairs = []
            if refs == -1 and refe == -1:
                # print(f"{qname}:{s}-{e} cannot be placed. Skipping..", file=sys.stderr)
                continue
            elif refs == -1:
                # print(f"Clipped start for {qname}:{s}-{e}. Skipping..", file=sys.stderr)
                for q, r in alpairs:
                    local_alpairs.append((q, r))
                    if r == refe:
                        break
            elif refe == -1:
                # print(f"Clipped end for {qname}:{s}-{e}. Skipping..", file=sys.stderr)
                for q, r in alpairs[::-1]:
                    local_alpairs.append((q, r))
                    if r == refs:
                        break
                local_alpairs.reverse()
            else:
                last_r = refs - 1
                for q, r in alpairs:
                    if r is None:
                        add_flag = last_r >= refs and last_r <= refe
                    else:
                        add_flag = refs <= r and r <= refe
                        last_r = r
                    if add_flag:
                        local_alpairs.append((q, r))

            if local_alpairs[0][1] == None or local_alpairs[-1][1] == None:
                # we have a clip
                clips_1 += 1
                continue

            # 3 we extract the 100 pairs preceding the region
            pre_alpairs = []
            add_flag = False
            n = 0
            for q, r in alpairs[::-1]:
                if add_flag and n < maxw:
                    pre_alpairs.append((q, r))
                    n += 1
                if r == local_alpairs[0][1]:
                    add_flag = True
            pre_alpairs.reverse()

            # 4 we extract the 100 pairs following the region
            post_alpairs = []
            add_flag = False
            n = 0
            for q, r in alpairs:
                if add_flag and n < maxw:
                    post_alpairs.append((q, r))
                    n += 1
                if r == local_alpairs[-1][1]:
                    add_flag = True

            # 5 we get the unique kmer in the upstream and downstream 100bp region
            (prek_q, prek_r), prekmerseq = get_unique_kmer(
                pre_alpairs, qseq, chroms[chrom], k, True)  # True for shorter
            (postk_q, postk_r), postkmerseq = get_unique_kmer(
                post_alpairs, qseq, chroms[chrom], k, False)  # False for shorter

            if prek_r == None:
                prek_r = local_alpairs[0][1]
            if postk_r == None:
                postk_r = local_alpairs[-1][1]

            # FIXME: sometimes the kmers are not correctly placed (None on read) - this should not happen: we have to extend over 100bp until a clean kmer
            if prek_q == None or postk_q == None:
                clips_2 += 1
                # print(f"{qname}:{s}-{e} cannot be extended on read. Skipping..", file=sys.stderr)
                continue
            if prek_r == None or postk_r == None:
                clips_2 += 1
                # print(f"{qname}:{s}-{e} cannot be extended on reference. Skipping..", file=sys.stderr)
                continue

            extsupersfss.add(SFS(chrom, qname, prek_r, postk_r + k,
                             prek_q, postk_q + k, qseq[prek_q:postk_q + k + 1]))
        n_al += 1
        if n_al % 1000 == 0:
            eprint(f"Parsed {n_al} alignments.")

    # NOTE: extsupersfss may contain multiple SFSs from the same read that cover the same variation. in any case, all of these will produce the same extended superstring

    eprint("Clustering..")
    tree = IntervalTree()
    for sfs in extsupersfss:
        overlaps = tree.overlap(sfs.rs, sfs.re+1)
        if len(overlaps) == 0:
            tree.addi(sfs.rs, sfs.re+1, [sfs])
        else:
            mins = sfs.rs
            maxe = sfs.re+1
            data = [sfs]
            for overlap in overlaps:
                mins = min(mins, overlap.begin)
                maxe = max(maxe, overlap.end)
                data.extend(overlap.data)
                tree.remove(overlap)
            tree.addi(mins, maxe, data)
    clusters = []
    for interval in tree:
        c = Cluster(interval.data[0].chrom)
        for sfs in interval.data:
            c.add(sfs)
        clusters.append(c)

    eprint(
        f"Analyzing {len(clusters)} clusters from {len(extsupersfss)} regions..")

    for cluster in clusters:
        if len(cluster) < minsupp:
            small_cl += 1
            continue

        chrom = cluster.chrom
        global_cluster = Cluster(cluster.chrom)
        reads = set()
        for sfs in cluster:
            reads.add(sfs.ridx)

        cov = 0
        for aln in bam.fetch(chrom, cluster.s, cluster.e):
            if aln.is_supplementary or aln.is_unmapped or aln.is_secondary:
                continue
            cov += 1
            qname = aln.query_name
            if qname not in reads:
                continue
            qseq = aln.query_sequence
            alpairs = aln.get_aligned_pairs()
            qs, qe = -1, -1
            for q, r in alpairs[::-1]:
                if q == None or r == None:
                    continue
                if r <= cluster.s:
                    qs = q
                    break
            for q, r in alpairs:
                if q == None or r == None:
                    continue
                if r >= cluster.e:
                    qe = q
                    break
            # TODO: get only prefix/suffix
            if qs == -1 or qe == -1:
                clips_3 += 1
                # print(f"Skipping {qname}:{qs}-{qe} due to clips", file=sys.stderr)
            else:
                global_cluster.add(
                    SFS(chrom, f"{qname}:{qs}-{qe}", cluster.s, cluster.e, qs, qe, qseq[qs:qe+1]))
        if len(global_cluster) >= minsupp:
            global_cluster.set_coverage(cov)
            global_cluster.dump()
            extcl += 1
        else:
            small_extcl += 1
    eprint(f"Dumped {extcl} clusters.")
    print("", file=sys.stderr)
    print(f"Clips #1 (SFS): {clips_1}", file=sys.stderr)
    print(f"Clips #2 (extension): {clips_1}", file=sys.stderr)
    print(f"Clips #3 (after extension): {clips_1}", file=sys.stderr)
    print(f"Small clusters: {small_cl}", file=sys.stderr)
    print(f"Small extended clusters: {small_extcl}", file=sys.stderr)


if __name__ == "__main__":
    main()
