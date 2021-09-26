#include <algorithm> 

#include "clipper.hpp"

using namespace std ;

Clipper::Clipper(const string& _chrom, const vector<Clip>& _clips) {
    chrom = _chrom ;
    clips = _clips ;
}

vector<Clip> Clipper::remove_duplicates(const vector<Clip> &clips) {
    vector<Clip> unique_clips;
    set<string> qnames;
    for (const Clip &clip : clips) {
        if (qnames.find(clip.name) == qnames.end()) {
            qnames.insert(clip.name);
            unique_clips.push_back(clip);
        }
    }
    return unique_clips;
}

vector<Clip> Clipper::combine(const vector<Clip> &clips) {
    vector<Clip> comb_clips;

    // we first cluster by breakpoints
    map<uint, vector<Clip>> clips_dict;
    for (const Clip &c : clips)
        clips_dict[c.p].push_back(c);
    // we then merge
    for (map<uint, vector<Clip>>::const_iterator it = clips_dict.begin(); it != clips_dict.end(); ++it) {
        uint max_l = 0;
        for (const Clip &c : it->second) {
            if (c.l > max_l) {
                max_l = c.l;
            }
        }
        comb_clips.push_back(Clip("", it->first, max_l, it->second.front().starting, it->second.size()));
    }
    return comb_clips;
}

vector<Clip> Clipper::filter_lowcovered(const vector<Clip> &clips, const uint w) {
    vector<Clip> filt_clips;
    for (const Clip &c : clips) {
        if (c.w >= w) {
            filt_clips.push_back(c);
        }
    }

    return filt_clips;
}

// Cluster clips by proximity
vector<Clip> Clipper::cluster(const vector<Clip> &clips, uint r) {
    vector<Clip> clusters;

    map<uint, Clip> clusters_by_pos;
    for (const Clip &c : clips) {
        bool found = false;
        for (map<uint, Clip>::iterator it = clusters_by_pos.begin(); it != clusters_by_pos.end(); ++it) {
            if (it->first - r <= c.p && c.p <= it->first + r) {
                found = true;
                it->second.l = max(it->second.l, c.l);
                it->second.w += c.w;
            }
        }

        if (!found) {
            clusters_by_pos[c.p] = c;
        }
    }

    for (map<uint, Clip>::iterator it = clusters_by_pos.begin(); it != clusters_by_pos.end(); ++it) {
        clusters.push_back(it->second);
    }
    return clusters;
}

vector<Clip> Clipper::filter_tooclose_clips(const vector<Clip> &clips, interval_tree_t<int> &vartree) {
    vector<Clip> fclips;

    for (const Clip &c : clips) {
        if (vartree.overlap_find({c.p, c.p + 1}) == std::end(vartree)) {
            fclips.push_back(c);
        }
    }

    return fclips;
}


void Clipper::call(const string &chrom_seq, interval_tree_t<int> &vartree) {
    vector<Clip> rclips;
    vector<Clip> lclips;
    for (const Clip &clip : clips) {
        if (clip.starting) {
            lclips.push_back(clip);
        } else {
            rclips.push_back(clip);
        }
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

    if (lclips.empty() || rclips.empty()) {
        return;
    }

    std::sort(lclips.begin(), lclips.end()) ;
    std::sort(rclips.begin(), rclips.end());

    for (const Clip &lc : lclips) {
        // we get the closest right clip
        // FIXME linear search is slow
        Clip rc;
        for (const Clip &rc_ : rclips) {
            if (lc.p < rc_.p) {
                rc = rc_;
                break;
            }
        }

        if (rc.w == 0) {
            continue;
        }

        if (abs((int)rc.p - (int)lc.p) < 1000) {
            // if a variation is in between the two clipped breakpoints
            // if vartree.overlaps(lS.s, rS.e):
            //     continue
            uint s = lc.w > rc.w ? lc.p : rc.p;
            uint l = max(lc.l, rc.l);
            string refbase = chrom_seq.substr(s, 1);
            uint w = max(lc.w, rc.w);

            svs.push_back(SV("INS", chrom, s, refbase, "<INS>", w, 0, 0, 0, true, l));
        }
    }

    for (const Clip &rc : rclips) {
        // we get the closest right clip
        Clip lc;
        for (const Clip &lc_ : lclips) {
            if (rc.p < lc_.p) {
                lc = lc_;
                break;
            }
        }

        if (lc.w == 0) {
            continue ;
        }

        if (lc.p - rc.p >= 2000 && lc.p - rc.p <= 50000) {
            uint s = rc.p;
            uint l = lc.p - rc.p + 1;
            string refbase = chrom_seq.substr(s, 1);
            uint w = max(lc.w, rc.w);
            if (w >= 5) {
                svs.push_back(SV("DEL", chrom, s, refbase, "<DEL>", w, 0, 0, 0, true, l));
            }
        }
    }
}
