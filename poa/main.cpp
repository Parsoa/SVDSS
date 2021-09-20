#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <regex>
#include <zlib.h>

#include "kseq.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "abpoa.h"
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#include "IntervalTree.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

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

struct Fragment
{
    int s;
    int e;
    string seq;

    Fragment(int s_, int e_, const string &seq_)
    {
        s = s_;
        e = e_;
        seq = seq_;
    }
};

struct Cluster
{
    string chrom;
    int s;
    int e;
    int cov;
    vector<string> seqs;

    Cluster(const string &chrom_, uint s_, uint e_, uint cov_)
    {
        chrom = chrom_;
        s = s_;
        e = e_;
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
};

struct SV
{
    string type;
    string chrom;
    int s;
    string refall;
    string altall;
    int e;
    int l;
    int w;
    int cov;
    int nv;
    int score;
    string gt;

    SV() { l = 0; }

    SV(const string type_, const string &chrom_, uint s_, const string &refall_, const string &altall_, const uint w_, const uint cov_, int l_, int nv_, int score_)
    {
        type = type_;
        chrom = chrom_;
        s = s_;
        refall = refall_;
        altall = altall_;
        e = s + refall.size() - 1;
        l = l_;
        w = w_;
        cov = cov_;
        score = score_;
        nv = nv_;
        genotype();
    }

    void genotype()
    {
        float p = (float)w / (float)cov;
        if (p<=0.1)
            gt = "0/0";
        else if (0.1 < p && p < 0.9)
            gt = "0/1";
        else
            gt = "1/1";
        // gt = './.';
    }
};

ostream &operator<<(ostream &os, const SV &sv)
{
    os << sv.chrom << "\t"
       << sv.s << "\t"
       << "."
       << "\t"
       << sv.refall << "\t"
       << sv.altall << "\t"
       << "."
       << "\t"
       << "PASS"
       << "\t"
       // INFO
       << "VARTYPE=SV;"
       << "SVTYPE=" << sv.type << ";"
       << "SVLEN=" << sv.l << ";"
       << "END=" << sv.e << ";"
       << "WEIGHT=" << sv.w << ";"
       << "COV=" << sv.cov << ";"
       << "NV=" << sv.nv << ";"
       << "ASCORE=" << sv.score
       << "\t"
       // -
       << "GT"
       << "\t"
       << sv.gt;
    return os;
}

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
    // o << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl;
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
        string info[4];
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

// def get_unique_kmer(alpairs, qseq, rseq, k, fromend=True):
//     kmers = {}
//     Is = list(range(0, len(alpairs) - k + 1))
//     for i in Is:
//         kmer = alpairs[i:i+k]
//         if any([q is None or r is None for (q, r) in kmer]):
//             # we want clean kmers only - ie placed kmers, no insertions or deletions
//             continue
//         # qkmer_seq = qseq[kmer[0][0]:kmer[-1][0]+1]
//         rkmer_seq = rseq[kmer[0][1]:kmer[-1][1]+1]
//         # note: after read reconstruction these kmers should be the same

//         kmers[rkmer_seq] = kmers[rkmer_seq] + 1 if rkmer_seq in kmers else 1

//     if fromend:
//         Is = Is[::-1]
//     last_kmer = [(None, None)]
//     last_rkmer_seq = ""
//     for i in Is:
//         kmer = alpairs[i:i+k]
//         if any([q is None or r is None for (q, r) in kmer]):
//             # we want clean kmers only - ie placed kmers, no insertions or deletions
//             continue
//         last_kmer = kmer
//         rkmer_seq = rseq[kmer[0][1]:kmer[-1][1]+1]
//         last_rkmer_seq = rkmer_seq
//         if kmers[rkmer_seq] == 1:
//             return kmer[0], rkmer_seq
//     return last_kmer[0], last_rkmer_seq

int main_extend(int argc, char *argv[])
{
    char *fa_path = argv[1];
    char *sfs_path = argv[2];
    char *bam_path = argv[3];

    // some hardcoded parameters FIXME
    uint k = 7;
    uint maxw = 100;
    uint minsupp = 2;

    // some stats on clips and clusters
    uint clips_1 = 0;
    uint clips_2 = 0;
    uint clips_3 = 0;
    uint cl = 0;
    uint small_cl = 0;
    uint small_extcl = 0;
    uint extcl = 0;

    // --- REFERENCE ---
    unordered_map<string, string> refs;
    gzFile fasta_in = gzopen(fa_path, "r");
    kseq_t *reference = kseq_init(fasta_in);
    int l;
    while ((l = kseq_read(reference)) >= 0)
        refs[reference->name.s] = reference->seq.s;
    kseq_destroy(reference);
    gzclose(fasta_in);

    // --- SFS ---
    map<string, vector<SFS>> SFSs = parse_sfsfile(sfs_path);

    // --- BAM ---
    // eprint("Parsing alignments..")
    // extsupersfss = set()
    samFile *bam = hts_open(bam_path, "r");
    hts_idx_t *bamidx = hts_idx_load(bam_path, HTS_FMT_BAI);
    bam_hdr_t *bamhdr = sam_hdr_read(bam);
    bam1_t *aln = bam_init1();
    uint n_al = 0;
    while (sam_read1(bam, bamhdr, aln) > 0)
    {
        if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY)
            continue;

        char *qname = bam_get_qname(aln);
        if (SFSs.find(qname) == SFSs.end())
            continue;

        char *chrom = bamhdr->target_name[aln->core.tid];
        uint qlen = aln->core.l_qseq;

        uint8_t *q = bam_get_seq(aln);
        char *qseq = (char *)malloc(qlen + 1);
        for (uint i = 0; i < qlen; ++i)
            qseq[i] = seq_nt16_str[bam_seqi(q, i)];
        qseq[qlen] = '\0';

        vector<pair<int, int>> alpairs = get_aligned_pairs(aln);

        for (SFS &sfs : SFSs.at(qname))
        {
            uint s = sfs.s;
            uint e = sfs.s + sfs.l - 1;
            // 1 we get the first bases on each direction (downstream and upstream the supersting) that have been correctly placed
            int refs = -1, refe = -1;
            for (int i = alpairs.size() - 1; i >= 0; --i)
            {
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (q == -1 || r == -1)
                    continue;
                if (q <= s)
                {
                    refs = r;
                    break;
                }
            }
            for (int i = 0; i < alpairs.size(); ++i)
            {
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (q == -1 || r == -1)
                    continue;
                if (q >= s)
                {
                    refe = r;
                    break;
                }
            }

            // 2 we extract the local alignment of the region of interest
            vector<pair<int, int>> local_alpairs;
            if (refs == -1 && refe == -1)
                continue;
            else if (refs == -1)
            {
                for (int i = 0; i < alpairs.size(); ++i)
                {
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    local_alpairs.push_back(make_pair(q, r));
                    if (r == refe)
                        break;
                }
            }
            else if (refe == -1)
            {
                for (int i = alpairs.size() - 1; i >= 0; --i)
                {
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    local_alpairs.push_back(make_pair(q, r));
                    if (r == refs)
                        break;
                }
                reverse(local_alpairs.begin(), local_alpairs.end());
            }
            else
            {
                int last_r = refs - 1;
                bool add_flag;
                for (int i = 0; i < alpairs.size(); ++i)
                {
                    int q = alpairs[i].first;
                    int r = alpairs[i].second;
                    if (r == -1)
                        add_flag = last_r >= refs && last_r <= refe;
                    else
                    {
                        add_flag = refs <= r && r <= refe;
                        last_r = r;
                    }
                    if (add_flag)
                        local_alpairs.push_back(make_pair(q, r));
                }
            }

            if (local_alpairs.front().first == -1 || local_alpairs.back().second == -1)
            {
                // we have a clip
                ++clips_1;
                continue;
            }

            // 3 we extract the 100 pairs preceding the region
            vector<pair<int, int>> pre_alpairs;
            bool add_flag = false;
            int n = 0;
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
                reverse(pre_alpairs.begin(), pre_alpairs.end());
            }

            // 4 we extract the 100 pairs following the region
            vector<pair<int, int>> post_alpairs;
            add_flag = false;
            n = 0;
            for (int i = 0; i < alpairs.size(); ++i)
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
            // 5 we get the unique kmer in the upstream and downstream 100bp region
            //     (prek_q, prek_r), prekmerseq = get_unique_kmer(
            //         pre_alpairs, qseq, chroms[chrom], k, True)  # True for shorter
            //     (postk_q, postk_r), postkmerseq = get_unique_kmer(
            //         post_alpairs, qseq, chroms[chrom], k, False)  # False for shorter

            //     if prek_r == None:
            //         prek_r = local_alpairs[0][1]
            //     if postk_r == None:
            //         postk_r = local_alpairs[-1][1]

            // FIXME: sometimes the kmers are not correctly placed (None on read) - this should not happen: we have to extend over 100bp until a clean kmer
            //     if prek_q == None or postk_q == None:
            //         clips_2 += 1
            //         # print(f"{qname}:{s}-{e} cannot be extended on read. Skipping..", file=sys.stderr)
            //         continue
            //     if prek_r == None or postk_r == None:
            //         clips_2 += 1
            //         # print(f"{qname}:{s}-{e} cannot be extended on reference. Skipping..", file=sys.stderr)
            //         continue

            //     extsupersfss.add(SFS(chrom, qname, prek_r, postk_r + k,
            //                      prek_q, postk_q + k, qseq[prek_q:postk_q + k + 1]))
            ++n_al;
            if (n_al % 1000 == 0)
                cerr << "Parsed " << n_al << " alignments." << endl;
        }
    }
}

int main(int argc, char *argv[])
{
    char *fa_path = argv[1];  // "ref.fa";
    char *bed_path = argv[2]; // "21.extended-regions.sorted.bed";
    string out_prefix = argv[3];

    uint minw = 2;
    uint mind = 15;

    ofstream osam(out_prefix + ".sam", ofstream::out);
    ofstream ovcf(out_prefix + ".vcf", ofstream::out);

    // --- REFERENCE
    unordered_map<string, string> refs;
    gzFile fasta_in = gzopen(fa_path, "r");
    kseq_t *reference = kseq_init(fasta_in);
    int l;
    while ((l = kseq_read(reference)) >= 0)
        refs[reference->name.s] = reference->seq.s;
    kseq_destroy(reference);
    gzclose(fasta_in);

    // --- EXTENDED REGIONS
    uint n_regions = 0;
    vector<Cluster> clusters;
    string line;
    ifstream inf(bed_path);
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

            clusters.push_back(Cluster(chrom, s, e, cov));
            for (i = 0; i < cl_size; ++i)
                clusters.back().add(seqs[i]);
            n_regions += cl_size;
        }
    }
    cerr << "Loaded " << clusters.size() << " clusters (" << n_regions << " extended supestrings)." << endl;

    // --- SAM header
    osam << "@HD\tVN:1.4" << endl;
    for (const auto &chrom : refs)
        osam << "@SQ\tSN:" << chrom.first << "\t"
             << "LN:" << chrom.second.size() << endl;

    vector<SV> all_svs;
    uint cl_idx = 0;
    for (const Cluster cl : clusters)
    {
        ++cl_idx;
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

        // --- Sorting clusters by #sequences
        sort(clusters_by_len.begin(), clusters_by_len.end(), [](const Cluster &c1, const Cluster &c2) -> bool
             { return c1.size() >= c2.size(); });

        // --- Analyze first 2 most weighted clusters
        for (uint i = 0; i < min((int)clusters_by_len.size(), 2); ++i)
        {
            Cluster c = clusters_by_len[i];
            if (c.size() < minw)
                continue;

            vector<SV> svs; // svs on current cluster
            string ref = refs[c.chrom].substr(c.s, c.e - c.s + 1);
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
                        svs.push_back(SV("INS", c.chrom, rpos, refs[c.chrom].substr(rpos - 1, 1), cons.substr(cpos, l), c.size(), c.cov, l, nv, score));
                    cpos += l;
                }
                else if (op == 'D')
                {
                    if (l > 25)
                        svs.push_back(SV("DEL", c.chrom, rpos, refs[c.chrom].substr(rpos - 1, l), refs[c.chrom].substr(rpos - 1, 1), c.size(), c.cov, -l, nv, score));
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
    // vector<Interval<uint, SV>> intervals;
    // for (const SV &sv : all_svs)
    //     intervals.push_back(Interval<uint, SV>(sv.s, sv.e+1, sv));
    // IntervalTree<uint, SV> tree = IntervalTree<uint, SV>(intervals);

    // -- VCF output
    print_vcf_header(refs, ovcf);
    for (const SV &sv : all_svs)
        // if (sv.gt.compare("0/0") != 0)
        ovcf << sv << endl;

    osam.close();
    ovcf.close();

    return 0;
}