#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <regex>
#include <zlib.h>

#include "kseq.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "abpoa.h"
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#include "interval_tree.hpp"

#include "sv.hpp"
#include "clipler.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace lib_interval_tree;

// AaCcGgTtNn ... ==> 0,1,2,3,4 ...
// BbDdEeFf   ... ==> 5,6,7,8 ...
unsigned char _char26_table[256] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 0, 5, 1, 6, 7, 8, 2, 9, 10, 11, 12, 13, 14, 4, 15,
    16, 17, 18, 19, 3, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26,
    26, 0, 5, 1, 6, 7, 8, 2, 9, 10, 11, 12, 13, 14, 4, 15,
    16, 17, 18, 19, 3, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26};

struct ExtSFS
{
    string chrom;
    string qname;
    int s;
    int e;
    // string seq;

    ExtSFS(const string &chrom_, const string &qname_, int s_, int e_ /*, const string &seq_*/)
    {
        chrom = chrom_;
        qname = qname_;
        s = s_;
        e = e_;
        // seq = seq_;
    }
};

struct Cluster
{
    string chrom;
    int s;
    int e;
    int cov;
    vector<string> seqs;

    Cluster(const string &chrom_, uint s_, uint e_, uint cov_ = 0)
    {
        chrom = chrom_;
        s = s_;
        e = e_;
        cov = cov_;
    }

    void set_cov(uint cov_)
    {
        cov = cov_;
    }

    void add(const string &seq)
    {
        seqs.push_back(seq);
    }

    int get_len() const
    {
        uint l = 0;
        uint n = 0;
        for (const string &seq : seqs)
        {
            ++n;
            l += seq.size();
        }
        return l / n;
    }

    vector<string> get_seqs() const
    {
        return seqs;
    }

    uint size() const
    {
        return seqs.size();
    }

    void dump(ofstream &o) const
    {
        o << chrom << " " << s << " " << e << " " << cov << " " << size();
        for (const string &seq : seqs)
            o << " " << seq;
        o << endl;
    }
};

void print_vcf_header(const unordered_map<string, string> &ref, ofstream &o)
{
    // TODO: improve and fix
    o << "##fileformat=VCFv4.2" << endl;
    // print("##fileDate=", date.today().strftime("%Y%m%d"), sep="")
    o << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz" << endl;
    for (unordered_map<string, string>::const_iterator it = ref.begin(); it != ref.end(); ++it)
        o << "##contig=<ID=" << it->first << ",length=" << it->second.size() << ">" << endl;
    o << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
    o << "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant class\">" << endl;
    o << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">" << endl;
    o << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
    o << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
    o << "##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number of alignments supporting this record\">" << endl;
    o << "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus\">" << endl;
    o << "##INFO=<ID=ASCORE,Number=1,Type=Integer,Description=\"Alignment score\">" << endl;
    o << "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of variations on same consensus\">" << endl;
    o << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl;
    o << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    o << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEFAULT" << endl;
}

struct SFS
{
    uint s;
    uint l;
    uint c;
    bool isreversed;

    SFS()
    {
        s = 0;
        l = 0;
        c = 0;
        isreversed = false;
    }

    SFS(uint s_, uint l_, uint c_, bool isreversed_)
    {
        s = s_;
        l = l_;
        c = c_;
        isreversed = isreversed_;
    }

    void reverse(uint p) { s = p - s - l; }
};

map<string, vector<SFS>> parse_sfsfile(const string &sfs_path)
{
    map<string, vector<SFS>> SFSs;
    string line;
    ifstream inf(sfs_path);
    if (inf.is_open())
    {
        string info[5];
        string read_name;
        while (getline(inf, line))
        {
            stringstream ssin(line);
            int i = 0;
            while (ssin.good() && i < 5)
            {
                ssin >> info[i++];
            }
            if (info[0].compare("*") != 0)
            {
                read_name = info[0];
                SFSs[read_name] = vector<SFS>();
            }
            //TODO
            SFSs[read_name].push_back(SFS(stoi(info[1]), stoi(info[2]), stoi(info[3]), true));
        }
    }
    return SFSs;
}

// Partial reimplementation of
// https://github.com/pysam-developers/pysam/blob/6ad0a57ef9c9b05d1492e10228ca7bccb5c7b30e/pysam/libcalignedsegment.pyx#L1867
// TODO run additional tests
vector<pair<int, int>> get_aligned_pairs(bam1_t *aln)
{
    vector<pair<int, int>> result;
    uint pos = aln->core.pos;
    uint qpos = 0;
    uint32_t *cigar = bam_get_cigar(aln);
    for (uint k = 0; k < aln->core.n_cigar; ++k)
    {
        uint32_t e = *(cigar + k);
        uint op = bam_cigar_op(e);   // e & 15
        uint l = bam_cigar_oplen(e); // e >> 4
        // char opc = bam_cigar_opchr(op);
        if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF)
        {
            for (uint i = pos; i < pos + l; ++i)
            {
                result.push_back(make_pair(qpos, i));
                ++qpos;
            }
            pos += l;
        }
        else if (op == BAM_CINS or op == BAM_CSOFT_CLIP)
        {
            for (uint i = pos; i < pos + l; ++i)
            {
                result.push_back(make_pair(qpos, -1));
                ++qpos;
            }
        }
        else if (op == BAM_CDEL)
        {
            for (uint i = pos; i < pos + l; ++i)
            {
                result.push_back(make_pair(-1, i));
            }
            pos += l;
        }
        else if (op == BAM_CHARD_CLIP)
        {
            // advances neither
        }
        else if (op == BAM_CREF_SKIP)
        {
            for (uint i = pos; i < pos + l; ++i)
            {
                result.push_back(make_pair(-1, i));
            }
            pos += l;
        }
        else if (op == BAM_CPAD)
        {
            // raise NotImplementedError(
            //     "Padding (BAM_CPAD, 6) is currently not supported. "
            //     "Please implement. Sorry about that.")
        }
    }
    return result;
}

pair<int, int> get_unique_kmer(const vector<pair<int, int>> &alpairs, const string &refseq, const uint k, const bool fromend)
{
    if (alpairs.size() < k)
        return make_pair(-1, -1);

    map<string, uint> kmers;
    string kmer_seq;
    for (uint i = 0; i < alpairs.size() - k + 1; ++i)
    {
        vector<pair<int, int>> kmer(alpairs.begin() + i, alpairs.begin() + i + k);
        bool skip = false;
        for (const pair<int, int> &p : kmer)
        {
            if (p.first == -1 || p.second == -1)
            {
                // we want clean kmers only - ie placed kmers, no insertions or deletions
                skip = true;
                break;
            }
        }
        if (skip)
            continue;
        kmer_seq = refseq.substr(kmer.front().second, k);
        // note: after read reconstruction, if the kmer is placed, than the kmer on the read and the kmer on the refernece should be the same
        kmers[kmer_seq] += 1;
    }

    vector<pair<int, int>> last_kmer;
    last_kmer.push_back(make_pair(-1, -1));
    for (uint i = 0; i < alpairs.size() - k + 1; ++i)
    {
        int new_i = i;
        if (fromend)
            new_i = alpairs.size() - k - i;
        vector<pair<int, int>> kmer(alpairs.begin() + new_i, alpairs.begin() + new_i + k);
        bool skip = false;
        for (const pair<int, int> &p : kmer)
        {
            if (p.first == -1 || p.second == -1)
            {
                // we want clean kmers only - ie placed kmers, no insertions or deletions
                skip = true;
                continue;
            }
        }
        if (skip)
            continue;
        last_kmer = kmer;
        kmer_seq = refseq.substr(kmer.front().second, k);
        if (kmers[kmer_seq] == 1)
            break;
    }
    return last_kmer.front();
}

int main_extend(char *fa_path, char *sfs_path, char *bam_path, char *out_prefix)
{
    string extsfss_path = string(out_prefix) + ".extsfss.bed";
    string clusters_path = string(out_prefix) + ".clusters.bed";
    string clips_path = string(out_prefix) + ".clips.bed";
    ofstream oextsfss(extsfss_path, ofstream::out);
    ofstream oclusters(clusters_path, ofstream::out);
    ofstream oclips(clips_path, ofstream::out);

    // some hardcoded parameters FIXME
    uint k = 7;
    uint maxw = 100;
    uint minsupp = 2;

    // --- REFERENCE ---
    unordered_map<string, string> reference;
    gzFile fasta_in = gzopen(fa_path, "r");
    kseq_t *ref = kseq_init(fasta_in);
    int l;
    while ((l = kseq_read(ref)) >= 0)
        reference[ref->name.s] = ref->seq.s;
    kseq_destroy(ref);
    gzclose(fasta_in);

    // --- SFS ---
    map<string, vector<SFS>> SFSs = parse_sfsfile(sfs_path);

    // --- BAM ---
    for (const auto &refseq : reference)
    {
        string chrom = refseq.first;
        cerr << "Processing chromosome " << chrom << ".. " << endl;

        // some stats
        uint skip_1 = 0;      // SFS skipped since no first/last base can be placed from read alignment (should be rare)
        uint skip_2 = 0;      // SFS skipped since it couldn't be extended
        uint skip_3 = 0;      // SFS skipped since reads starts/ends inside a cluster
        uint small_cl = 0;    // number of cluster (before clustering) with low support
        uint extcl = 0;       // number of extended clusters (after clustering)
        uint small_extcl = 0; // number of extended clusters (after clustering) with low support

        vector<ExtSFS> extsfss;
        samFile *bam = hts_open(bam_path, "r");
        hts_idx_t *bamidx = hts_idx_load(bam_path, HTS_FMT_BAI);
        bam_hdr_t *bamhdr = sam_hdr_read(bam);
        bam1_t *aln = bam_init1();
        hts_itr_t *itr = sam_itr_querys(bamidx, bamhdr, chrom.c_str());
        uint n_al = 0;
        while (sam_itr_next(bam, itr, aln) > 0)
        {
            if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY)
                continue;

            char *qname = bam_get_qname(aln);
            if (SFSs.find(qname) == SFSs.end())
                continue;

            uint aln_start = aln->core.pos;
            uint aln_end = bam_endpos(aln);
            uint32_t *cigar = bam_get_cigar(aln);

            vector<pair<int, int>> alpairs = get_aligned_pairs(aln);

            // NOTE: we may have more sfs on a clipped read, but all of them will produce the same clip
            pair<uint, uint> lclip;
            pair<uint, uint> rclip;
            for (SFS &sfs : SFSs.at(qname))
            {
                int s = sfs.s;
                int e = sfs.s + sfs.l - 1;

                // 1 we get the first bases on each direction (downstream and upstream the supersting) that have been correctly placed
                int refs = -1, refe = -1;
                for (int i = alpairs.size() - 1; i >= 0; --i)
                {
                    // checking from end to rightmost placed base
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (q == -1 || r == -1)
                        continue;
                    if (q < s) // <= seems more correct to me but using < we are more flexible
                    {
                        refs = r;
                        break;
                    }
                }
                for (uint i = 0; i < alpairs.size(); ++i)
                {
                    // checking from begin to leftmost placed base
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (q == -1 || r == -1)
                        continue;
                    if (q > e) // >= seems more correct but using > we are more flexible
                    {
                        refe = r;
                        break;
                    }
                }

                // 2 we extract the local alignment of the region of interest
                vector<pair<int, int>> local_alpairs;
                if (refs == -1 && refe == -1)
                {
                    // we couldn't place the first and the last base, so we skip this - otherwise we'll end up considering the entire read
                    ++skip_1;
                    continue;
                }
                else if (refs == -1)
                {
                    // // we couldn't place the first base, so we fill with all the starting pairs we have (ie we have a prefix of the alignment)..
                    // for (uint i = 0; i < alpairs.size(); ++i)
                    // {
                    //     int q = alpairs[i].first;
                    //     int r = alpairs[i].second;
                    //     local_alpairs.push_back(make_pair(q, r));
                    //     if (r == refe)
                    //         break;
                    // }
                    // ..actually we just extract the clips from the alignment
                    uint op = bam_cigar_op(*(cigar + 0));
                    uint l = bam_cigar_oplen(*(cigar + 0));
                    if (op == BAM_CSOFT_CLIP)
                        lclip = make_pair(aln_start, l);
                    // we skip this SFS
                    continue;
                }
                else if (refe == -1)
                {
                    // // we couldn't place the last base, so we fill with all the ending pairs we have (ie we have a suffix of the alignment)..
                    // for (int i = alpairs.size() - 1; i >= 0; --i)
                    // {
                    //     int q = alpairs[i].first;
                    //     int r = alpairs[i].second;
                    //     local_alpairs.push_back(make_pair(q, r));
                    //     if (r == refs)
                    //         break;
                    // }
                    // reverse(local_alpairs.begin(), local_alpairs.end());
                    // ..actually we just extract the clips from the alignment
                    uint op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
                    uint l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
                    if (op == BAM_CSOFT_CLIP)
                        rclip = make_pair(aln_end, l);
                    // we skip this SFS
                    continue;
                }
                else
                {
                    // we placed the first and last base, so we extract the alignment (ie a substring)
                    int last_r = refs - 1;
                    bool add_flag;
                    for (uint i = 0; i < alpairs.size(); ++i)
                    {
                        int q = alpairs[i].first;
                        int r = alpairs[i].second;
                        if (r == -1)
                            add_flag = refs <= last_r && last_r <= refe;
                        else
                        {
                            add_flag = refs <= r && r <= refe;
                            last_r = r;
                        }
                        if (add_flag)
                            local_alpairs.push_back(make_pair(q, r));

                        // We break when we found a placed base at or after the reference end
                        if (q != -1 && r != -1 && r >= refe)
                            break;
                    }
                }

                // only if we have been able to place the SFS..
                // 3 ..we extract the maxw (100) pairs preceding the region of interest
                vector<pair<int, int>> pre_alpairs;
                bool add_flag = false;
                uint n = 0;
                for (int i = alpairs.size() - 1; i >= 0; --i)
                {
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (add_flag && n < maxw)
                    {
                        pre_alpairs.push_back(make_pair(q, r));
                        ++n;
                    }
                    if (r == local_alpairs.front().second)
                    {
                        add_flag = true;
                    }
                }
                reverse(pre_alpairs.begin(), pre_alpairs.end());

                // 4 we extract the maxw (100) pairs following the region of interest
                vector<pair<int, int>> post_alpairs;
                add_flag = false;
                n = 0;
                for (uint i = 0; i < alpairs.size(); ++i)
                {
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (add_flag and n < maxw)
                    {
                        post_alpairs.push_back(make_pair(q, r));
                        ++n;
                    }
                    if (r == local_alpairs.back().second)
                        add_flag = true;
                }

                // 5 we get the unique kmer in the upstream and downstream maxw-bp regions
                pair<int, int> prekmer = get_unique_kmer(pre_alpairs, reference[chrom], k, true);    // true for first kmer found (shorter cluster)
                pair<int, int> postkmer = get_unique_kmer(post_alpairs, reference[chrom], k, false); // false for first kmer found (shorter cluster)

                // if we couldn't place a kmer, we just get the entire region
                if (prekmer.first == -1 || prekmer.second == -1)
                {
                    prekmer.first = local_alpairs.front().first;
                    prekmer.second = local_alpairs.front().second;
                }
                if (postkmer.first == -1 || postkmer.second == -1)
                {
                    postkmer.first = local_alpairs.back().first;
                    postkmer.second = local_alpairs.back().second;
                }

                // if also the entire region is not correctly placed, then we skip it
                // NOTE: I think we can solve this by increasing maxw
                if (prekmer.first == -1 || prekmer.second == -1 || postkmer.first == -1 || postkmer.second == -1)
                {
                    ++skip_2;
                    continue;
                }
                // otherwise, we store it as an extended SFS
                // FIXME: understand why this is happening (chr16 on full giab genome)
                if ((uint)prekmer.second > postkmer.second + k)
                    cerr << "Error on " << qname << ". SFS starting at " << sfs.s << " (length " << sfs.l << ")." << endl;
                else
                    extsfss.push_back(ExtSFS(string(chrom), string(qname), prekmer.second, postkmer.second + k));
            }
            // we dump clips on the current read
            // TODO: store total coverage
            if (lclip.second > 0)
                oclips << chrom << "\t" << lclip.first << "\t" << lclip.first + 1 << "\t" << lclip.second << "\t" << qname << "\t"
                       << "L" << endl;
            if (rclip.second > 0)
                oclips << chrom << "\t" << rclip.first << "\t" << rclip.first + 1 << "\t" << rclip.second << "\t" << qname << "\t"
                       << "R" << endl;
            ++n_al;
            if (n_al % 10000 == 0)
                cerr << "Parsed " << n_al << " alignments." << endl;
        }

        // NOTE: extsupersfss may contain multiple SFSs from the same read that cover the same variation. in any case, all of these will produce the same extended superstring

        // --- CLUSTERING
        interval_tree_t<int> tree;
        for (const ExtSFS &sfs : extsfss)
        {
            oextsfss << sfs.chrom << "\t" << sfs.s << " " << sfs.e << " " << sfs.qname << endl;
            vector<pair<int, int>> overlaps;
            tree.overlap_find_all({sfs.s, sfs.e}, [&overlaps](auto iter)
                                  {
                                      overlaps.push_back(make_pair(iter->low(), iter->high()));
                                      return true;
                                  });
            if (overlaps.empty())
            {
                tree.insert({sfs.s, sfs.e});
            }
            else
            {
                int mins = sfs.s;
                int maxe = sfs.e;
                for (const pair<int, int> overlap : overlaps)
                {
                    mins = min(mins, overlap.first);
                    maxe = max(maxe, overlap.second);
                }
                tree.insert({mins, maxe});
            }
        }
        tree.deoverlap();

        map<pair<int, int>, vector<ExtSFS>> clusters;
        for (const ExtSFS &sfs : extsfss)
        {
            auto overlap = tree.overlap_find({sfs.s, sfs.e});
            clusters[make_pair(overlap->low(), overlap->high())].push_back(sfs);
        }

        // --- CLUSTERS CLEANING
        cerr << "Analyzing " << clusters.size() << " clusters from " << extsfss.size() << " extSFSs.." << endl;

        for (const auto &cluster : clusters)
        {
            string chrom = cluster.second.front().chrom;
            set<string> reads;
            int cluster_s = numeric_limits<int>::max();
            int cluster_e = 0;
            for (const ExtSFS &esfs : cluster.second)
            {
                cluster_s = min(cluster_s, esfs.s);
                cluster_e = max(cluster_e, esfs.e);
                reads.insert(esfs.qname);
            }
            uint cluster_size = reads.size();
            if (cluster_size < minsupp)
            {
                ++small_cl;
                continue;
            }

            Cluster global_cluster = Cluster(chrom, cluster_s, cluster_e);

            uint cov = 0;
            aln = bam_init1();
            string region = chrom + ":" + to_string(cluster_s) + "-" + to_string(cluster_e);
            hts_itr_t *itr = sam_itr_querys(bamidx, bamhdr, region.c_str());
            while (sam_itr_next(bam, itr, aln) > 0)
            {
                if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY)
                    continue;
                ++cov; // FIXME: this cov takes into account also reads starting or ending inside the cluster (maybe we should skip those?)

                char *qname = bam_get_qname(aln);
                if (reads.find(qname) == reads.end())
                    continue;

                uint qlen = aln->core.l_qseq;
                uint8_t *q = bam_get_seq(aln);
                char *qseq = (char *)malloc(qlen + 1);
                for (uint i = 0; i < qlen; ++i)
                    qseq[i] = seq_nt16_str[bam_seqi(q, i)];
                qseq[qlen] = '\0';

                vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
                int qs = -1, qe = -1;
                // getting starting and ending positions on read sequence aligning to start/end of cluster
                for (int i = alpairs.size() - 1; i >= 0; --i)
                {
                    // finding starting position
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (q == -1 || r == -1)
                        continue;
                    if (r <= cluster_s)
                    {
                        qs = q;
                        break;
                    }
                }
                for (uint i = 0; i < alpairs.size(); ++i)
                {
                    // finding ending position
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (q == -1 || r == -1)
                        continue;
                    if (r >= cluster_e)
                    {
                        qe = q;
                        break;
                    }
                }

                if (qs == -1 || qe == -1)
                    // reads starts or ends inside the cluster
                    // TODO: get only remaining prefix/suffix? but this may make POA and realignment harder
                    ++skip_3;
                else
                {
                    string seq(qseq, qs, qe - qs + 1);
                    global_cluster.add(seq);
                }
            }
            if (global_cluster.size() >= minsupp)
            {
                global_cluster.set_cov(cov);
                global_cluster.dump(oclusters);
                ++extcl;
            }
            else
            {
                ++small_extcl;
            }
        }
        cerr << "Dumped " << extcl << " clusters." << endl
             << endl;
        cerr << "Skipped due to unplaced alignment: " << skip_1 << endl;
        cerr << "Skipped due to unplaced extension: " << skip_2 << endl;
        cerr << "Skipped due to start/end inside extended cluster: " << skip_3 << endl;
        cerr << "Skipped extended clusters (due to support): " << small_extcl << endl;
    }

    oextsfss.close();
    oclips.close();
    oclusters.close();

    return 0;
}

int main_call(char *fa_path, char *out_prefix)
{
    string clusters_path = string(out_prefix) + ".clusters.bed"; // FIXME: hardcoded as in main_extend
    string clips_path = string(out_prefix) + ".clips.bed";       // FIXME: hardcoded as in main_extend
    string sam_path = string(out_prefix) + ".consensus.sam";
    string vcf_path = string(out_prefix) + ".variations.vcf";

    ofstream osam(sam_path, ofstream::out);
    ofstream ovcf(vcf_path, ofstream::out);

    // some hardcoded parameters FIXME
    uint minw = 2;
    uint mind = 15;

    // --- REFERENCE
    unordered_map<string, string> reference;
    gzFile fasta_in = gzopen(fa_path, "r");
    kseq_t *refseq = kseq_init(fasta_in);
    int l;
    while ((l = kseq_read(refseq)) >= 0)
        reference[refseq->name.s] = refseq->seq.s;
    kseq_destroy(refseq);
    gzclose(fasta_in);

    // --- EXTENDED REGIONS
    uint n_regions = 0;
    map<string, vector<Cluster>> clusters;
    string line;
    ifstream inf(clusters_path);
    if (inf.is_open())
    {
        string line;
        while (getline(inf, line))
        {
            stringstream ssin(line);

            string chrom, token;
            ssin >> chrom;
            ssin >> token;
            uint s = stoi(token);
            ssin >> token;
            uint e = stoi(token);
            ssin >> token;
            uint cov = stoi(token);
            ssin >> token;
            uint cl_size = stoi(token);

            string seqs[cl_size + 1];
            uint i = 0;
            while (ssin.good() && i < cl_size + 1)
                ssin >> seqs[i++];

            clusters[chrom].push_back(Cluster(chrom, s, e, cov));
            for (i = 0; i < cl_size; ++i)
                clusters[chrom].back().add(seqs[i]);
            n_regions += cl_size;
        }
    }
    inf.close();
    cerr << "Loaded " << clusters.size() << " clusters (" << n_regions << " extSFSs)." << endl;

    // --- SAM header
    osam << "@HD\tVN:1.4" << endl;
    for (const auto &chrom : reference)
        osam << "@SQ\tSN:" << chrom.first << "\t"
             << "LN:" << chrom.second.size() << endl;

    // --- VCF header
    print_vcf_header(reference, ovcf);

    for (const auto &chr_clusters : clusters)
    {
        string chrom = chr_clusters.first;
        cerr << "Calling on chromosome " << chrom << ".." << endl;
        vector<SV> all_svs;
        for (const Cluster &cl : chr_clusters.second)
        {
            if (cl.size() < minw)
                continue;

            // --- Splitting cluster by sequence length
            vector<Cluster> clusters_by_len;
            for (const string &seq : cl.get_seqs())
            {
                uint i;
                for (i = 0; i < clusters_by_len.size(); ++i)
                {
                    if (abs((int)clusters_by_len[i].get_len() - (int)seq.size()) <= mind)
                        break;
                }
                if (i == clusters_by_len.size())
                    clusters_by_len.push_back(Cluster(cl.chrom, cl.s, cl.e, cl.cov));
                clusters_by_len[i].add(seq);
            }

            // --- Sorting clusters by #sequences to get first 2 most weighted clusters
            int i_max1 = -1;
            int i_max2 = -1;
            uint v_max1 = 0;
            uint v_max2 = 0;
            for (uint i = 0; i < clusters_by_len.size(); ++i)
            {
                if (clusters_by_len[i].size() > v_max1)
                {
                    v_max2 = v_max1;
                    i_max2 = i_max1;
                    v_max1 = clusters_by_len[i].size();
                    i_max1 = i;
                }
                else if (clusters_by_len[i].size() > v_max2)
                {
                    v_max2 = clusters_by_len[i].size();
                    i_max2 = i;
                }
            }
            vector<int> maxs ({i_max1, i_max2});

            // NOTE: I tried the following but it crashed on chr8 (big cluster 43822762-43837766)
            // sort(clusters_by_len.begin(), clusters_by_len.end(), [](Cluster &c1, Cluster &c2) -> bool
            //      { return c1.size() > c2.size(); });
            // for (int i = 0; i < min((int)clusters_by_len.size(), 2); ++i)

            for (const int i : maxs)
            {
                if (i == -1)
                    continue;
                Cluster c = clusters_by_len[i];
                if (c.size() < minw)
                    continue;

                vector<SV> svs; // svs on current cluster
                string ref = reference[c.chrom].substr(c.s, c.e - c.s + 1);
                vector<string> seqs = c.get_seqs();

                // --- POA
                uint n_seqs = seqs.size();
                abpoa_t *ab = abpoa_init();
                abpoa_para_t *abpt = abpoa_init_para();
                abpt->disable_seeding = 1;
                abpt->align_mode = 0; // global
                abpt->out_msa = 0;
                abpt->out_cons = 1;
                abpt->out_gfa = 0;
                // abpt->w = 6, abpt->k = 9;
                // abpt->min_w = 10; // minimizer-based seeding and partition
                abpt->is_diploid = 1;
                // abpt->min_freq = 0.3;
                abpt->progressive_poa = 0;
                abpoa_post_set_para(abpt);

                int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
                uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
                for (uint i = 0; i < n_seqs; ++i)
                {
                    seq_lens[i] = seqs[i].size();
                    bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
                    for (int j = 0; j < seq_lens[i]; ++j)
                        bseqs[i][j] = _char26_table[(int)seqs[i][j]];
                    // cerr << seqs[i] << endl;
                }
                uint8_t **cons_seq;
                int **cons_cov, *cons_l, cons_n = 0;
                uint8_t **msa_seq;
                int msa_l = 0;
                abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);

                if (cons_l[0] == 0)
                    continue;
                string cons = "";
                for (int j = 0; j < cons_l[0]; ++j)
                    cons += "ACGTN-"[cons_seq[0][j]];
                // cerr << ">ref" << "\n" << ref << endl;
                // cerr << ">cons" << "\n" << cons << endl;

                // --- Local realignment
                parasail_result_t *result = NULL;
                result = parasail_nw_trace_striped_16(cons.c_str(), cons.size(), ref.c_str(), ref.size(), 10, 1, &parasail_blosum62);
                parasail_cigar_t *cigar = parasail_result_get_cigar(result, cons.c_str(), cons.size(), ref.c_str(), ref.size(), NULL);
                string cigar_str = parasail_cigar_decode(cigar);
                int score = result->score;

                // Dumping SAM
                osam << c.chrom << ":" << c.s + 1 << "-" << c.e + 1 << "_" << c.size() << "\t"
                     << "0"
                     << "\t" << c.chrom << "\t" << c.s + 1 << "\t"
                     << "60"
                     << "\t" << cigar_str << "\t"
                     << "*"
                     << "\t"
                     << "0"
                     << "\t"
                     << "0"
                     << "\t" << cons << "\t"
                     << "*" << endl;

                // -- Parsing CIGAR
                vector<pair<uint, char>> cigar_pairs;
                regex r("([0-9]+)([MIDNSHPX=])");
                regex_iterator<string::iterator> rit(cigar_str.begin(), cigar_str.end(), r);
                regex_iterator<string::iterator> rend;
                uint nv = 0;
                while (rit != rend)
                {
                    int l = stoi(rit->str().substr(0, rit->str().size() - 1));
                    char op = rit->str().substr(rit->str().size() - 1, 1)[0];
                    cigar_pairs.push_back(make_pair(l, op));
                    if (l > 25 && (op == 'I' || op == 'D'))
                        ++nv;
                    ++rit;
                }

                // -- Extracting SVs
                uint rpos = c.s; // position on reference
                uint cpos = 0;   // position on consensus
                for (const auto cigar_pair : cigar_pairs)
                {
                    uint l = cigar_pair.first;
                    char op = cigar_pair.second;
                    if (op == '=' || op == 'M')
                    {
                        rpos += l;
                        cpos += l;
                    }
                    else if (op == 'I')
                    {
                        if (l > 25)
                            svs.push_back(SV("INS", c.chrom, rpos, reference[c.chrom].substr(rpos - 1, 1), cons.substr(cpos, l), c.size(), c.cov, l, nv, score));
                        cpos += l;
                    }
                    else if (op == 'D')
                    {
                        if (l > 25)
                            svs.push_back(SV("DEL", c.chrom, rpos, reference[c.chrom].substr(rpos - 1, l), reference[c.chrom].substr(rpos - 1, 1), c.size(), c.cov, -l, nv, score));
                        rpos += l;
                    }
                }
                parasail_cigar_free(cigar);
                parasail_result_free(result);

                // --- combine SVs on same consensus ---
                // TODO: add a flag to do this
                vector<SV> merged_svs;
                // TODO improve this
                SV ins;
                SV del;
                for (const auto &sv : svs)
                {
                    if (sv.l < 0)
                    {
                        if (del.l == 0)
                        {
                            del = sv;
                            continue;
                        }
                        else
                        {
                            del.l += sv.l;
                            del.refall += sv.refall;
                            del.e += (-sv.l);
                        }
                    }
                    else
                    {
                        if (ins.l == 0)
                        {
                            ins = sv;
                            continue;
                        }
                        else
                        {
                            ins.l += sv.l;
                            ins.altall += sv.altall;
                        }
                    }
                }
                if (ins.l != 0)
                    merged_svs.push_back(ins);
                if (del.l != 0)
                    merged_svs.push_back(del);

                // -- Combine svs with same length (maybe useless now - only if diploid mode)
                vector<SV> comb_svs;
                for (const SV &msv : merged_svs)
                {
                    bool newsv_flag = true;
                    for (SV &csv : comb_svs)
                    {
                        if (abs(csv.l - msv.l) <= 10)
                        {
                            csv.w += msv.w;
                            newsv_flag = false;
                        }
                    }
                    if (newsv_flag)
                    {
                        comb_svs.push_back(msv);
                    }
                }
                for (const SV &sv : comb_svs)
                    all_svs.push_back(sv);

                // -- ABPOA CLEANING --
                if (cons_n)
                {
                    for (int i = 0; i < cons_n; ++i)
                    {
                        free(cons_seq[i]);
                        // free(cons_cov[i]);
                    }
                    free(cons_seq);
                    // free(cons_cov);
                    free(cons_l);
                }
                if (msa_l)
                {
                    for (uint i = 0; i < n_seqs; ++i)
                        free(msa_seq[i]);
                    free(msa_seq);
                }
                for (uint i = 0; i < n_seqs; ++i)
                    free(bseqs[i]);
                free(bseqs);
                free(seq_lens);
                abpoa_free(ab);
                abpoa_free_para(abpt);
                // -- END ABPOA CLEANING --
            }
        }

        // --- SV clustering
        // this is currently done in the cluster_sv.py python script

        // -- VCF output
        interval_tree_t<int> sv_tree;
        for (const SV &sv : all_svs)
        {
            // if (sv.gt.compare("0/0") != 0)
            ovcf << sv << endl;
            sv_tree.insert({sv.s, sv.e});
        }

        // -- CLIPS
        // Clipler clipler(chrom, clips_path);
        // clipler.call(reference[chrom], sv_tree);
        // for (const SV &sv : clipler.osvs)
        //     ovcf << sv << endl;
    }

    osam.close();
    ovcf.close();

    return 0;
}

int main(int argc, char *argv[])
{
    char *fa_path = argv[1];
    char *sfs_path = argv[2];
    char *bam_path = argv[3];
    char *out_prefix = argv[4];

    // TODO
    // * merge and polish the two methods
    // * keep clusters and clips also in memory
    main_extend(fa_path, sfs_path, bam_path, out_prefix);
    main_call(fa_path, out_prefix);
}