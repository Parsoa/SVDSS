#include "insdeller.hpp"

Insdeller::Insdeller(const string &chrom_, samFile *sfs_bam_, bam_hdr_t *sfs_bamhdr_, hts_idx_t *sfs_bamindex_, samFile *read_bam_, bam_hdr_t *read_bamhdr_, hts_idx_t *read_bamindex_)
{
    chrom = chrom_;
    sfs_bam = sfs_bam_;
    sfs_bamhdr = sfs_bamhdr_;
    sfs_bamindex = sfs_bamindex_;

    read_bam = read_bam_;
    read_bamhdr = read_bamhdr_;
    read_bamindex = read_bamindex_;
}

list<Cluster> Insdeller::pcluster()
{
    bam1_t *aln = bam_init1();
    hts_itr_t *itr = sam_itr_querys(sfs_bamindex, sfs_bamhdr, chrom.c_str());

    list<Cluster> clusters; // this will store the clusters
    clusters.push_back(Cluster(chrom));

    while (sam_itr_next(sfs_bam, itr, aln) > 0)
    {
        string sfs_name(bam_get_qname(aln));
        string qname = sfs_name.substr(0, sfs_name.find(".")); // keep just read name (no positions on read)
        uint rs = aln->core.pos;
        uint re = bam_endpos(aln);

        bool has_i = false;
        bool has_d = false;
        uint32_t *cigar = bam_get_cigar(aln);
        for (uint i = 0; i < aln->core.n_cigar; ++i)
        {
            uint op = bam_cigar_op(*(cigar + i));

            has_i |= op == BAM_CINS;
            has_d |= op == BAM_CDEL;
        }

        if (!has_i && !has_d)
            continue;
        uint ft = 0;
        if (has_i && has_d)
            ft = 2;
        else if (has_i)
            ft = 1;

        Cluster &c = clusters.back();
        if (c.empty())
        {
            c.add_fragment(Fragment(qname, rs, re, ft));
        }
        else
        {
            if (c.back().e + 500 > rs) // overlapping
            {
                c.add_fragment(Fragment(qname, rs, re, ft));
            }
            else
            {
                clusters.push_back(Cluster(chrom));
                clusters.back().add_fragment(Fragment(qname, rs, re, ft));
            }
        }
    }

    return clusters;
}

// Type clustering - 2 clusters with all fragments with only I and only D or 1 cluster with both (if we have at least one mixed fragment)
list<Cluster> Insdeller::tcluster(const Cluster &c)
{
    list<Cluster> tc;
    bool has_both = false;
    for (const Fragment &f : c)
        has_both |= f.t == 2;

    if (has_both)
        tc.push_back(c);
    else
    {
        tc.push_back(Cluster(chrom));
        tc.push_back(Cluster(chrom));
        for (const Fragment &f : c)
        {
            Cluster &c = f.t == 0 ? tc.front() : tc.back();
            c.add_fragment(f);
        }
    }
    return tc;
}

Cluster Insdeller::extend(const Cluster &c)
{
    Cluster ext_c(chrom);

    int rs = c.s; // reference start
    int re = c.e; // reference end

    set<string> reads_in_cluster; // these are the reads we are interested in
    for (const Fragment f : c)
        reads_in_cluster.insert(f.name);

    string region = chrom + ":" + to_string(rs) + "-" + to_string(re);
    hts_itr_t *itr = sam_itr_querys(read_bamindex, read_bamhdr, region.c_str());
    bam1_t *aln = bam_init1();
    while (sam_itr_next(read_bam, itr, aln) > 0)
    {
        string qname = bam_get_qname(aln);
        if (reads_in_cluster.find(qname) == reads_in_cluster.end())
            continue;
        int qs = -1; // read start
        int qe = -1; // read end
        vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
        for (uint i = 0; i < alpairs.size(); ++i)
        {
            int q = alpairs[i].first;
            int r = alpairs[i].second;
            if (r == rs)
            {
                qs = q;
            }
            else if (r == re)
            {
                qe = q;
            }
        }

        // --- this should never occur since we already re-aligned the superstrings (but I could be wrong)
        if (qs == -1)
            qs = alpairs.front().first;
        if (qe == -1)
            qe = alpairs.back().first;
        if (qs == -1 || qe == -1)
        {
            cerr << "[W] no reference base on " << chrom << ":" << rs << "-" << re << " (" << qname << ")" << endl;
            continue;
        }
        // ---

        // add flanking
        uint flank = 10; // FIXME: hardcoded
        qs -= flank;
        qe += flank;

        // FIXME: may over/underflow
        qs = max(0, qs);

        uint qlen = aln->core.l_qseq;
        uint8_t *q = bam_get_seq(aln);
        char *qseq = (char *)malloc(qlen + 1); // FIXME: avoid this malloc everytime
        for (uint i = 0; i < qlen; ++i)
            qseq[i] = seq_nt16_str[bam_seqi(q, i)];
        qseq[qlen] = '\0';
        string subseq(qseq, qs, qe - qs + 1);
        Fragment f(qname, rs - flank, re + flank, 0, subseq); // FIXME 0 is good? // CHECKME the size of the cluster is the initial size - even if we split the cluster
        ext_c.add_fragment(f);
    }
    return ext_c;
}

list<Cluster> Insdeller::scluster(const Cluster &cluster) // CHECKME
{
    list<Cluster> clusters;

    // greedy choice: longer fragment is representative #1
    Fragment repr1;
    for (const Fragment &f : cluster)
        if (f.size() > repr1.size())
            repr1 = f;

    // the other representative will be the fragment with the lowest similarity with repr1
    Fragment repr2;
    double min_d = 100.0;
    for (const Fragment &f : cluster)
    {
        double d = rapidfuzz::fuzz::ratio(repr1.seq, f.seq);
        if (d < 60 && d < min_d) // CHECKME 60 is hardcoded
        {
            min_d = d;
            repr2 = f;
        }
    }

    if (repr2.size() == 0)
    {
        clusters.push_back(cluster);
    }
    else
    {
        clusters.push_back(Cluster(cluster.chrom));
        clusters.push_back(Cluster(cluster.chrom));
        double d1, d2;
        for (const Fragment &f : cluster)
        {
            d1 = rapidfuzz::fuzz::ratio(repr1.seq, f.seq);
            d2 = rapidfuzz::fuzz::ratio(repr2.seq, f.seq);
            if (d1 <= d2)
                clusters.front().add_fragment(f);
            else
                clusters.back().add_fragment(f);
        }
    }
    return clusters;
}

CIGAR Insdeller::align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs;
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    ts = (uint8_t *)malloc(tl);
    qs = (uint8_t *)malloc(ql);
    for (i = 0; i < tl; ++i)
        ts[i] = _nt4_table[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i)
        qs[i] = _nt4_table[(uint8_t)qseq[i]];
    // ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);

    ksw_extd(0, ql, qs, tl, ts, 5, mat, 50, 1, 50, 0, -1, -1, 0, &ez);

    vector<pair<uint, char>> cigar(ez.n_cigar);
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
        cigar[i] = make_pair(ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
    int score = ez.score;
    free(ez.cigar);
    free(ts);
    free(qs);
    return CIGAR(cigar, score);
}

list<SV> Insdeller::extract(const Cluster &c, const string &ref, ofstream &o)
{
    string reference(ref, c.s, c.e - c.s + 1);
    string consensus = c.poa();
    CIGAR cigar = align(reference.c_str(), consensus.c_str(), 1, -3, 100, 0); // TODO: adjust scores/penalties. I want less gaps
    cigar.fixclips();

    // SAM output
    o << chrom << ":" << c.s + 1 << "-" << c.e + 1 << ":" << c.size() << "\t"
      << 0 << "\t"
      << c.chrom << "\t"
      << c.s + 1 << "\t"
      << 60 << "\t"
      << cigar.to_str() << "\t"
      << "*"
      << "\t"
      << 0 << "\t"
      << 0 << "\t"
      << consensus << "\t"
      << "*"
      << endl;
    // << "NM:i:" << 0 << endl;
    // << "AS:i:" << cigar.score << endl;

    // extracting insertions/deletions breakpoints
    uint rs = c.s + 1; // position on reference
    uint lrs = 0;      // position on local reference
    uint qs = 0;       // position on consensus
    string refall = "";
    string altall = "";

    list<SV> svs;
    for (uint i = 0; i < cigar.size(); ++i)
    {
        SV sv;
        switch (cigar[i].second)
        {
        case 'S':
            qs += cigar[i].first;
            break;
        case 'I':
            refall = reference.substr(lrs - 1, 1);
            altall = refall + consensus.substr(qs, cigar[i].first);
            sv = SV('I', c.chrom, rs - 1, refall, altall, c.size(), c.full_cov, cigar.ngaps, cigar.score);
            qs += cigar[i].first;
            break;
        case 'D':
            refall = reference.substr(lrs - 1, cigar[i].first + 1);
            altall = refall.substr(0, 1);
            sv = SV('D', c.chrom, rs - 1, refall, altall, c.size(), c.full_cov, cigar.ngaps, cigar.score);
            rs += cigar[i].first;
            lrs += cigar[i].first;
            break;
        case 'M':
            lrs += cigar[i].first;
            rs += cigar[i].first;
            qs += cigar[i].first;
            break;
        default:
            cerr << "Unknown CIGAR op " << cigar[i].second << endl;
            exit(1);
        }
        if (abs(sv.l) > 0)
            svs.push_back(sv);
    }
    return svs;
}

list<SV> Insdeller::merge_svs(const list<SV> &svs, const string &chrom_seq)
{
    vector<list<SV>> tsvs(2);
    list<SV> msvs;
    for (const auto &sv : svs)
    {
        uint i = sv.type.compare("DEL") == 0 ? 0 : 1;
        tsvs[i].push_back(sv);
    }
    // Note: we are sure SVs (of same type) are sorted by reference position
    for (uint i = 0; i < 2; ++i)
    {
        if (tsvs[i].empty())
            continue;
        string stype = tsvs[i].front().type;
        char type;
        string chrom = tsvs[i].front().chrom;
        uint s = tsvs[i].front().s;
        string refall = "";
        string altall = "";
        uint w = tsvs[i].front().w;
        uint cov = tsvs[i].front().cov;
        int ngaps = tsvs[i].front().ngaps;
        int score = tsvs[i].front().score;
        if (stype.compare("INS") == 0)
        {
            type = 'I';
            refall = tsvs[i].front().refall;
            for (const SV &sv : tsvs[i])
                altall += sv.altall;
        }
        else
        {
            type = 'D';
            altall = tsvs[i].front().altall;
            uint l = 0;
            for (const auto sv : tsvs[i])
                l += sv.refall.size();
            refall = chrom_seq.substr(s, l);
        }

        msvs.push_back(SV(type, chrom, s, refall, altall, w, cov, ngaps, score));
    }
    return msvs;
}

list<SV> Insdeller::dedup_svs(const list<SV> &svs)
{
    set<string> uidxs;
    list<SV> usvs;
    for (const SV &sv : svs)
    {
        if (uidxs.find(sv.idx) == uidxs.end())
        {
            uidxs.insert(sv.idx);
            usvs.push_back(sv);
        }
    }
    return usvs;
}

void Insdeller::call(const string &chrom_seq, ofstream &osam)
{
    list<Cluster> pclusters = pcluster();

    for (const Cluster &pc : pclusters)
    {
        if (pc.size() < 2)
            continue;

        list<Cluster> tclusters = tcluster(pc);
        for (const Cluster &tc : tclusters)
        {
            if (tc.size() < 2)
                continue;

            list<SV> cl_svs;

            Cluster extended_tc = extend(tc);
            list<Cluster> sclusters = scluster(extended_tc);
            for (Cluster &sc : sclusters)
            {
                if (sc.size() < 2)
                    continue;

                // get total number of reads mapping to this locus
                bam1_t *aln = bam_init1();
                string region = sc.chrom + ":" + to_string(sc.s) + "-" + to_string(sc.e + 1);
                hts_itr_t *itr = sam_itr_querys(read_bamindex, read_bamhdr, region.c_str());
                uint c = 0;
                while (sam_itr_next(read_bam, itr, aln) > 0)
                {
                    if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY)
                        continue;
                    ++c;
                }
                sc.set_full_coverage(c);

                list<SV> svs = extract(sc, chrom_seq, osam);
                list<SV> msvs = merge_svs(svs, chrom_seq);

                cl_svs.insert(cl_svs.end(), msvs.begin(), msvs.end());
            }

            if (!cl_svs.empty())
            {
                // CHECKME maybe useless
                list<SV> usvs = dedup_svs(cl_svs);

                for (const SV &sv : usvs)
                {
                    if (abs(sv.l) >= 30 && (sv.ngaps <= 2 || (sv.ngaps > 2 && sv.w > 10)))
                    {
                        osvs.push_back(sv);
                        vartree.insert({sv.s, sv.e});
                    }
                }
            }
        }
    }
}