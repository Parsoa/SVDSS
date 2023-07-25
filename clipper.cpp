#include "clipper.hpp"

Clipper::Clipper(const vector<Clip> &_clips) { clips = _clips; }

// vector<Clip> Clipper::remove_duplicates(const vector<Clip> &clips) {
//   vector<Clip> unique_clips;
//   unordered_map<string, int> qnames;
//   for (const Clip &clip : clips) {
//     if (qnames.find(clip.name) == qnames.end()) {
//       qnames[clip.name] = 0;
//       unique_clips.push_back(clip);
//     }
//   }
//   return unique_clips;
// }

// vector<Clip> Clipper::combine(const vector<Clip> &clips) {
//   int threads = 4;
//   vector<vector<Clip>> _p_combined_clips;
//   _p_combined_clips.resize(threads);
//   // we first cluster by breakpoints
//   unordered_map<string, unordered_map<uint, vector<Clip>>> clips_dict;
//   for (const Clip &c : clips) {
//     clips_dict[c.chrom][c.p].push_back(c);
//   }
// // we then merge
// #pragma omp parallel for num_threads(threads) schedule(static, 1)
//   for (int i = 0; i < chromosomes.size(); i++) {
//     int t = i % threads;
//     const string &chrom = chromosomes[i];
//     for (auto it = clips_dict[chrom].begin(); it != clips_dict[chrom].end();
//          ++it) {
//       uint max_l = 0;
//       for (const Clip &c : it->second) {
//         if (c.l > max_l) {
//           max_l = c.l;
//         }
//       }
//       Clip clip = Clip("", chrom, it->first, max_l,
//       it->second.front().starting,
//                        it->second.size());
//       _p_combined_clips[t].push_back(clip);
//     }
//   }
//   vector<Clip> combined_clips;
//   for (int i = 0; i < threads; i++) {
//     combined_clips.insert(combined_clips.begin(),
//     _p_combined_clips[i].begin(),
//                           _p_combined_clips[i].end());
//   }
//   return combined_clips;
// }

// vector<Clip> Clipper::filter_lowcovered(const vector<Clip> &clips,
//                                         const uint w) {
//   vector<Clip> filtered_clips;
//   for (const Clip &c : clips) {
//     if (c.w >= w) {
//       filtered_clips.push_back(c);
//     }
//   }
//   return filtered_clips;
// }

// // Cluster clips by proximity
// // TODO: this might be too slow
// vector<Clip> Clipper::cluster(const vector<Clip> &clips, uint r) {
//   vector<Clip> clusters;
//   map<uint, Clip> clusters_by_pos;
//   for (const Clip &c : clips) {
//     bool found = false;
//     for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
//          it != clusters_by_pos.end(); ++it) {
//       if (it->first - r <= c.p && c.p <= it->first + r) {
//         found = true;
//         it->second.l = max(it->second.l, c.l);
//         it->second.w += c.w;
//       }
//     }
//     if (!found) {
//       clusters_by_pos[c.p] = c;
//     }
//   }

//   for (map<uint, Clip>::iterator it = clusters_by_pos.begin();
//        it != clusters_by_pos.end(); ++it) {
//     clusters.push_back(it->second);
//   }
//   return clusters;
// }

// vector<Clip> Clipper::filter_tooclose_clips(const vector<Clip> &clips,
//                                             interval_tree_t<int> &vartree) {
//   vector<Clip> fclips;
//   for (const Clip &c : clips) {
//     if (vartree.overlap_find({c.p, c.p + 1}) == end(vartree)) {
//       fclips.push_back(c);
//     }
//   }
//   return fclips;
// }

// // find smallest right that is larger than query
// int binary_search(const vector<Clip> &clips, int begin, int end,
//                   const Clip &query) {
//   // for (int i = 0; i < clips.size(); i++) {
//   //     if (query.p < clips[i].p) {
//   //         return i ;
//   //     }
//   // }
//   // return -1 ;
//   if (begin > end || begin >= clips.size()) {
//     return -1;
//   }
//   int m = (begin + end) / 2;
//   if (clips[m].p == query.p) {
//     if (m + 1 < clips.size()) {
//       return m + 1;
//     } else {
//       return m;
//     }
//   } else if (clips[m].p > query.p) {
//     if (m > 0 && clips[m - 1].p < query.p) {
//       return m;
//     }
//     return binary_search(clips, begin, m - 1, query);
//   } else {
//     return binary_search(clips, m + 1, end, query);
//   }
// }

// void Clipper::call(int threads, interval_tree_t<int> &vartree) {
//   // lprint({"Predicting SVS from", to_string(clips.size()), "clipped SFS
//   on",
//   //         to_string(threads), "threads.."});
//   vector<Clip> rclips;
//   vector<Clip> lclips;
//   for (const Clip &clip : clips) {
//     if (clip.starting) {
//       lclips.push_back(clip);
//     } else {
//       rclips.push_back(clip);
//     }
//   }
//   // lprint({to_string(lclips.size()), "left clips."});
//   // lprint({to_string(rclips.size()), "right clips."});
//   // lprint({"Preprocessing clipped SFS.."});
// #pragma omp parallel for num_threads(2) schedule(static, 1)
//   for (int i = 0; i < 2; i++) {
//     if (i == 0) {
//       rclips = remove_duplicates(rclips);
//       rclips = combine(rclips);
//       rclips = filter_lowcovered(rclips, 2); // FIXME: hardcoded
//       rclips = filter_tooclose_clips(rclips, vartree);
//       rclips = cluster(rclips, 1000); // FIXME: hardcoded
//       sort(rclips.begin(), rclips.end());
//     } else {
//       lclips = remove_duplicates(lclips);
//       lclips = combine(lclips);
//       lclips = filter_lowcovered(lclips, 2); // FIXME: hardcoded
//       lclips = filter_tooclose_clips(lclips, vartree);
//       lclips = cluster(lclips, 1000); // FIXME: hardcoded
//       sort(lclips.begin(), lclips.end());
//     }
//   }
//   // lprint({to_string(lclips.size()), "left clips."});
//   // lprint({to_string(rclips.size()), "right clips."});
//   _p_svs.resize(threads);
//   if (lclips.empty() || rclips.empty()) {
//     return;
//   }
//   // lprint({"Predicting insertions.."});
// #pragma omp parallel for num_threads(threads) schedule(static, 1)
//   for (int i = 0; i < lclips.size(); i++) {
//     const Clip &lc = lclips[i];
//     int t = omp_get_thread_num();
//     string chrom = lc.chrom;
//     // we get the closest right clip
//     int r = binary_search(rclips, 0, rclips.size() - 1, lc);
//     if (r == -1) {
//       continue;
//     }
//     auto rc = rclips[r];
//     if (rc.w == 0) {
//       continue;
//     }

//     if (abs((int)rc.p - (int)lc.p) < 1000) {
//       uint s = lc.w > rc.w ? lc.p : rc.p;
//       uint l = max(lc.l, rc.l);
//       string refbase(chromosome_seqs[chrom] + s, 1);
//       uint w = max(lc.w, rc.w);
//       _p_svs[t].push_back(
//           SV("INS", chrom, s, refbase, "<INS>", w, 0, 0, 0, true, l));
//     }
//   }
//   // lprint({"Predicting deletions.."});
// #pragma omp parallel for num_threads(threads) schedule(static, 1)
//   for (int i = 0; i < rclips.size(); i++) {
//     const Clip &rc = rclips[i];
//     int t = omp_get_thread_num();
//     string chrom = rc.chrom;
//     // we get the closest right clip
//     int l = binary_search(lclips, 0, lclips.size() - 1, rc);
//     if (l == -1) {
//       continue;
//     }
//     auto lc = lclips[l];
//     if (lc.w == 0) {
//       continue;
//     }

//     if (lc.p - rc.p >= 2000 && lc.p - rc.p <= 50000) {
//       uint s = rc.p;
//       uint l = lc.p - rc.p + 1;
//       string refbase(chromosome_seqs[chrom] + s, 1);
//       uint w = max(lc.w, rc.w);
//       if (w >= 5) {
//         _p_svs[t].push_back(
//             SV("DEL", chrom, s, refbase, "<DEL>", w, 0, 0, 0, true, l));
//       }
//     }
//   }
// }
