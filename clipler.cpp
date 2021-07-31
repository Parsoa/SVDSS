#include "clipler.hpp"

list<Clip> Clipler::remove_duplicates(const list<Clip> &clips)
{
    list<Clip> unique_clips;
    set<string> qnames;
    for (const Clip &clip : clips)
    {
        if (qnames.find(clip.name) == qnames.end())
        {
            qnames.insert(clip.name);
            unique_clips.push_back(clip);
        }
    }
    return unique_clips;
}

list<Clip> Clipler::combine(const list<Clip> &clips)
{
    list<Clip> comb_clips;

    // we first cluster by breakpoints
    map<uint, list<Clip>> clips_dict;
    for (const Clip &c : clips)
        clips_dict[c.p].push_back(c);
    // we then merge
    for (map<uint, list<Clip>>::const_iterator it = clips_dict.begin(); it != clips_dict.end(); ++it)
    {
        uint max_l = 0;
        for (const Clip &c : it->second)
        {
            if (c.l > max_l)
            {
                max_l = c.l;
            }
        }
        comb_clips.push_back(Clip("", it->first, max_l, it->second.front().starting, it->second.size()));
    }
    return comb_clips;
}

list<Clip> Clipler::filter_lowcovered(const list<Clip> &clips, const uint w)
{
    list<Clip> filt_clips;
    for (const Clip &c : clips)
    {
        if (c.w >= w)
            filt_clips.push_back(c);
    }

    return filt_clips;
}

// Cluster clips by proximity
list<Clip> Clipler::cluster(const list<Clip> &clips, uint r)
{
    list<Clip> clusters;

    map<uint, Clip> clusters_by_pos;
    for (const Clip &c : clips)
    {
        bool found = false;
        for (map<uint, Clip>::iterator it = clusters_by_pos.begin(); it != clusters_by_pos.end(); ++it)
        {
            if (it->first - r <= c.p && c.p <= it->first + r)
            {
                found = true;
                it->second.l = max(it->second.l, c.l);
                it->second.w += c.w;
            }
        }

        if (!found)
        {
            clusters_by_pos[c.p] = c;
        }
    }

    for (map<uint, Clip>::iterator it = clusters_by_pos.begin(); it != clusters_by_pos.end(); ++it)
        clusters.push_back(it->second);
    return clusters;
}

list<Clip> Clipler::filter_tooclose_clips(const list<Clip> &clips, interval_tree_t<int> &vartree)
{
    list<Clip> fclips;

    for (const Clip &c : clips)
    {
        if (vartree.overlap_find({c.p, c.p + 1}) == std::end(vartree))
        {
            fclips.push_back(c);
        }
    }

    return fclips;
}

Clipler::Clipler(const string &chrom_, samFile *sfs_bam_, bam_hdr_t *sfs_bamhdr_, hts_idx_t *sfs_bamindex_)
{
    chrom = chrom_;
    sfs_bam = sfs_bam_;
    sfs_bamhdr = sfs_bamhdr_;
    sfs_bamindex = sfs_bamindex_;
}

list<Clip> Clipler::extract_clips()
{
    bam1_t *aln = bam_init1();
    hts_itr_t *itr = sam_itr_querys(sfs_bamindex, sfs_bamhdr, chrom.c_str());
    // string region = "1:102568879-102568949";
    // hts_itr_t *itr = sam_itr_querys(sfs_bamindex, sfs_bamhdr, region.c_str());

    list<Clip> clips;
    while (sam_itr_next(sfs_bam, itr, aln) > 0)
    {
        string sfs_name(bam_get_qname(aln));
        string qname = sfs_name.substr(0, sfs_name.find(".")); // keep just read name (no positions on read)
        uint rs = aln->core.pos;
        uint re = bam_endpos(aln);
        uint32_t *cigar = bam_get_cigar(aln);

        // first
        uint op = bam_cigar_op(*(cigar + 0));
        uint l = bam_cigar_oplen(*(cigar + 0));
        if (op == BAM_CSOFT_CLIP)
            clips.push_back(Clip(qname, rs, l, true));

        // last
        op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
        l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
        if (op == BAM_CSOFT_CLIP)
            clips.push_back(Clip(qname, re, l, false));
    }
    return clips;
}

void Clipler::call(const string &chrom_seq, interval_tree_t<int> &vartree)
{
    list<Clip> clips = extract_clips();

    list<Clip> rclips;
    list<Clip> lclips;
    for (const Clip &clip : clips)
    {
        if (clip.starting)
            lclips.push_back(clip);
        else
            rclips.push_back(clip);
    }

    rclips = remove_duplicates(rclips);
    lclips = remove_duplicates(lclips);
    lclips = combine(lclips);
    rclips = combine(rclips);
    lclips = filter_lowcovered(lclips, 2); // FIXME: hardcoded
    rclips = filter_lowcovered(rclips, 2); // FIXME: hardcoded
    lclips = filter_tooclose_clips(lclips, vartree);
    rclips = filter_tooclose_clips(rclips, vartree);
    lclips = cluster(lclips, 1000); // FIXME: hardcoded
    rclips = cluster(rclips, 1000); // FIXME: hardcoded

    if (lclips.empty() || rclips.empty())
        return;

    lclips.sort();
    rclips.sort();

    for (const Clip &lc : lclips)
    {
        // we get the closest right clip
        Clip rc;
        for (const Clip &rc_ : rclips)
        {
            if (lc.p < rc_.p)
            {
                rc = rc_;
                break;
            }
        }

        if (rc.w == 0)
            continue;

        if (abs((int)rc.p - (int)lc.p) < 1000) // FIXME: hardcoded
        {
            // if a variation is in between the two clipped breakpoints
            // if vartree.overlaps(lS.s, rS.e):
            //     continue
            uint s = lc.w > rc.w ? lc.p : rc.p;
            uint l = max(lc.l, rc.l);
            string refbase = chrom_seq.substr(s, 1);
            uint w = max(lc.w, rc.w);

            // TODO: get coverage of locus from bam
            osvs.push_back(SV('I', chrom, s, refbase, "<INS>", w, 0, 0, 0, true, l));
        }
    }

    for (const Clip &rc : rclips)
    {
        // we get the closest right clip
        Clip lc;
        for (const Clip &lc_ : lclips)
        {
            if (rc.p < lc_.p)
            {
                lc = lc_;
                break;
            }
        }

        if (lc.w == 0)
            continue;

        if (lc.p - rc.p >= 2000 && lc.p - rc.p <= 50000) // FIXME: hardcoded
        {
            uint s = rc.p;
            uint l = lc.p - rc.p + 1;
            string refbase = chrom_seq.substr(s, 1);
            uint w = max(lc.w, rc.w);

            // TODO: get coverage of locus from bam
            if (w >= 5)
                osvs.push_back(SV('D', chrom, s, refbase, "<DEL>", w, 0, 0, 0, true, l));
        }
    }
}