#include <ctime>

#include "extender.hpp"

using namespace std;
using namespace lib_interval_tree;

Extender::Extender(unordered_map<string, vector<SFS>> *_SFSs) {
  SFSs = _SFSs;
  config = Configuration::getInstance();
  min_w = config->min_cluster_weight;
}

void Extender::run(int _threads) {
  threads = _threads;
  batch_size = int(10000 / threads) * threads;
  bam_file = hts_open(config->bam.c_str(), "r");
  bam_index = sam_index_load(bam_file, config->bam.c_str());
  bam_header = sam_hdr_read(bam_file);
  bgzf_mt(bam_file->fp.bgzf, 8, 1);
  // extend reads
  _p_clips.resize(threads);
  _p_extended_sfs.resize(threads);
  extend_parallel();
  for (int i = 0; i < threads; i++) {
    // extended_sfs.insert(extended_sfs.begin(), _p_extended_sfs[i].begin(),
    // _p_extended_sfs[i].end()) ;
    for (const auto &extsfs : _p_extended_sfs[i]) {
      extended_sfs.push_back(extsfs);
    }
    clips.insert(clips.begin(), _p_clips[i].begin(), _p_clips[i].end());
  }
  lprint({"Extension complete.", to_string(clips.size()), " clipped SFS."});
  // build a separate interval tree for each chromosome
  cluster_no_interval_tree();
  // put all clusters in a single vector
  lprint({"Flattening interval clusters.."});
  map<pair<int, int>, ExtCluster> _ext_clusters;
  for (int i = 0; i < threads; i++) {
    for (const auto &cluster : _p_sfs_clusters[i]) {
      _ext_clusters.insert(
          make_pair(cluster.first, ExtCluster(cluster.second)));
    }
  }

  for (const auto &cluster : _ext_clusters)
    ext_clusters.push_back(cluster.second);

  // process each cluster separately
  _p_clusters.resize(threads);
  extract_sfs_sequences();
  // flatten clusters into a single vector
  for (int i = 0; i < threads; i++) {
    clusters.insert(clusters.begin(), _p_clusters[i].begin(),
                    _p_clusters[i].end());
  }

  ofstream clofile;
  clofile.open(config->workdir + "/clusters.txt");
  for (const auto &cluster : clusters) {
    // CHECKME: is ending position (cluster.first.second) inclusive? I'm
    // assuming yes
    clofile << cluster.chrom << ":" << cluster.s << "-" << cluster.e << "\t"
            << cluster.size();
    for (int i = 0; i < cluster.size(); ++i)
      clofile << "\t" << cluster.get_name(i) << ":" << cluster.get_seq(i);
    clofile << endl;
  }
  clofile.close();

  // merge POA alignments
  _p_svs.resize(threads);
  _p_alignments.resize(threads);
  call();
  for (int i = 0; i < threads; i++) {
    svs.insert(svs.begin(), _p_svs[i].begin(), _p_svs[i].end());
    alignments.insert(alignments.begin(), _p_alignments[i].begin(),
                      _p_alignments[i].end());
  }
  lprint({"Extracted", to_string(svs.size()), "SVs."});
  sam_close(bam_file);
  filter_sv_chains();
  // lprint({"Error stats: ", to_string(skip_1), " can't be mapped, ",
  // to_string(skip_2), " can't be extended ", to_string(skip_3), " anomalies.",
  // to_string(skip_4), " clipped."}) ;
}

pair<int, int> Extender::get_unique_kmers(const vector<pair<int, int>> &alpairs,
                                          const uint k, const bool from_end,
                                          string chrom) {
  if (alpairs.size() < k) {
    return make_pair(-1, -1);
  }

  map<string, int> kmers;
  string kmer_seq;
  int i = 0;
  while (i < alpairs.size() - k + 1) {
    bool skip = false;
    for (int j = i; j < i + k; j++) {
      if (alpairs[j].first == -1 || alpairs[j].second == -1) {
        // we want clean kmers only - ie placed kmers, no insertions or
        // deletions
        skip = true;
        i = j + 1; // jump to next possible start position
        break;
      }
    }
    if (skip) {
      continue;
    }
    string kmer(chromosome_seqs[chrom] + alpairs[i].second, k);
    kmers[kmer] += 1;
    i++;
  }
  pair<int, int> last_kmer = make_pair(-1, -1);
  i = 0;
  while (i < alpairs.size() - k + 1) {
    int offset = i;
    if (from_end) {
      offset = alpairs.size() - k - i;
    }
    bool skip = false;
    for (int j = offset; j < offset + k; j++) {
      if (alpairs[j].first == -1 || alpairs[j].second == -1) {
        skip = true;
        i += (j - offset); // jump to next possible start position
        break;
      }
    }
    if (skip) {
      i++;
      continue;
    }
    last_kmer = alpairs[offset];
    string kmer(chromosome_seqs[chrom] + alpairs[offset].second, k);
    if (kmers[kmer] == 1) {
      break;
    }
    i++;
  }
  return last_kmer;
}

// Parallelize within each chromosome
void Extender::extend_parallel() {
  lprint({"Extending superstrings on", to_string(threads), "threads.."});
  int p = 0;
  int b = 0;
  for (int i = 0; i < 2; i++) {
    bam_entries.push_back(vector<vector<bam1_t *>>(threads));
    for (int j = 0; j < threads; j++) {
      for (int k = 0; k < batch_size / threads; k++) {
        bam_entries[i][j].push_back(bam_init1());
      }
    }
  }
  load_batch_bam(threads, batch_size, p);
  time_t t;
  time(&t);
  bool should_load = true;
  bool should_process = true;
  bool should_terminate = false;
  bool loaded_last_batch = false;
  uint64_t u = 0;
  int num_reads = 0;
  bool printed = false;
#pragma omp parallel num_threads(config->threads + 1)
  while (!loaded_last_batch) {
// lprint({"Beginning batch", to_string(b + 1)});
#pragma omp single
    {
      for (int i = 0; i < threads; i++) {
        u += bam_entries[p][i].size();
      }
    }
#pragma omp for
    for (int i = 0; i < threads + 1; i++) {
      int t = omp_get_thread_num();
      if (t == 0) {
        // load next batch of entries
        if (should_load) {
          loaded_last_batch = !load_batch_bam(threads, batch_size, (p + 1) % 2);
        }
      } else {
        process_batch(bam_entries[p][t - 1], t - 1);
      }
    }
#pragma omp single
    {
      p += 1;
      p %= 2;
      b += 1;
      time_t s;
      time(&s);
      if (s - t == 0) {
        s += 1;
      }
      cerr << "[I] Processed batch " << std::left << b
           << ". Alignments so far: " << std::right << u
           << ". Alignments per second: " << u / (s - t)
           << ". Time: " << std::fixed << s - t << "\r";
      printed = true;
    }
  }
  if (printed) {
    cerr << endl;
  }
  // cleanup
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < threads; j++) {
      for (int k = 0; k < batch_size / threads; k++) {
        bam_destroy1(bam_entries[i][j][k]);
      }
    }
  }
  lprint({"Done."});
}

bool Extender::load_batch_bam(int threads, int batch_size, int p) {
  int i = 0;
  int n = 0;
  while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
    bam1_t *aln = bam_entries[p][n % threads][i];
    if (aln == nullptr) {
      break;
    }
    if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY ||
        aln->core.flag & BAM_FSECONDARY) {
      continue;
    }
    char *qname = bam_get_qname(aln);
    if (SFSs->find(qname) == SFSs->end()) {
      continue;
    }
    n += 1;
    if (n % threads == 0) {
      i += 1;
    }
    if (n == batch_size) {
      return true;
    }
  }
  if (n != 0 && n != batch_size) {
    for (int j = n % threads; j < threads; j++)
      bam_entries[p][j][i] = nullptr;
    for (int j = 0; j < n % threads; j++)
      if (i + 1 < bam_entries[p][j].size())
        bam_entries[p][j][i + 1] = nullptr;
  }
  // lprint({"Loaded", to_string(n), "BAM reads.."});
  return n != 0 ? true : false;
}

void Extender::process_batch(vector<bam1_t *> bam_entries, int index) {
  bam1_t *aln;
  for (int b = 0; b < bam_entries.size(); b++) {
    aln = bam_entries[b];
    if (aln == nullptr) {
      break;
    }
    extend_alignment(aln, index);
  }
}

void Extender::extend_alignment(bam1_t *aln, int index) {
  char *qname = bam_get_qname(aln);
  uint32_t *cigar = bam_get_cigar(aln);
  vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
  string chrom(bam_header->target_name[aln->core.tid]);
  if (chromosome_seqs.find(chrom) == chromosome_seqs.end()) {
    return;
  }

  // NOTE: we may have more sfs on a clipped read, but all of them will produce
  // the same clip
  pair<uint, uint> lclip;
  pair<uint, uint> rclip;
  int last_pos = 0;
  vector<ExtSFS> local_extended_sfs;
  for (const SFS &sfs : SFSs->at(qname)) {
    int s = sfs.s;
    int e = sfs.s + sfs.l - 1;
    int aln_start = -1;
    int aln_end = -1;
    vector<pair<int, int>> local_alpairs;
    // find start and end of SFS ibn read's alignment
    int refs = -1;
    int refe = -1;
    for (int i = last_pos; i < alpairs.size(); i++) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      if (q == -1 || r == -1) {
        continue;
      } else if (q < s) { // <= seems more correct to me but using < we are more
                          // flexible
        last_pos = i;
        refs = r;
        aln_start = i;
      } else if (q >
                 e) { // >= seems more correct but using > we are more flexible
        refe = r;
        aln_end = i;
        break;
      }
    }
    // We extract the local alignment of the region of interest
    if (refs == -1 && refe == -1) {
      // we couldn't place the first and the last base, so we skip this -
      // otherwise we'll end up considering the entire read
      ++skip_1;
      continue;
    } else if (refs == -1) {
      uint op = bam_cigar_op(*(cigar + 0));
      uint l = bam_cigar_oplen(*(cigar + 0));
      if (op == BAM_CSOFT_CLIP) {
        lclip = make_pair(aln->core.pos, l);
        skip_4++;
      }
      // we skip this SFS
      continue;
    } else if (refe == -1) {
      uint op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
      uint l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
      if (op == BAM_CSOFT_CLIP) {
        rclip = make_pair(bam_endpos(aln), l);
        skip_4++;
      }
      // we skip this SFS
      continue;
    } else {
      // we placed the first and last base, so we extract the alignment (ie a
      // substring)
      int last_r = refs - 1;
      for (uint i = aln_start; i <= aln_end; i++) {
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (r == -1) {
          if (refs <= last_r && last_r <= refe) {
            local_alpairs.push_back(make_pair(q, r));
          }
        } else {
          last_r = r;
          if (refs <= r && r <= refe) {
            local_alpairs.push_back(make_pair(q, r));
          }
        }
        // We break when we found a placed base at or after the reference end
        if (q != -1 && r != -1 && r >= refe) {
          break;
        }
      }
    }

    // only if we have been able to place the SFS..
    // 3 ..we extract the maxw (100) pairs preceding the region of interest
    vector<pair<int, int>> pre_alpairs;
    uint n = 0;
    for (int i = aln_start - 1; i >= 0; --i) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      if (n < maxw) {
        pre_alpairs.push_back(make_pair(q, r));
        n++;
      }
    }
    reverse(pre_alpairs.begin(), pre_alpairs.end());
    // 4 we extract the maxw (100) pairs following the region of interest
    vector<pair<int, int>> post_alpairs;
    n = 0;
    for (uint i = aln_end + 1; i < alpairs.size(); i++) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      if (n < maxw) {
        post_alpairs.push_back(make_pair(q, r));
        n++;
      }
    }

    // 5 we get the unique kmer in the upstream and downstream maxw-bp regions
    pair<int, int> prekmer =
        get_unique_kmers(pre_alpairs, kmer_size, true,
                         chrom); // true for first kmer found (shorter cluster)
    pair<int, int> postkmer =
        get_unique_kmers(post_alpairs, kmer_size, false,
                         chrom); // false for first kmer found (shorter cluster)

    // if we couldn't place a kmer, we just get the entire region
    if (prekmer.first == -1 || prekmer.second == -1) {
      prekmer.first = local_alpairs.front().first;
      prekmer.second = local_alpairs.front().second;
    }
    if (postkmer.first == -1 || postkmer.second == -1) {
      postkmer.first = local_alpairs.back().first;
      postkmer.second = local_alpairs.back().second;
    }
    // if also the entire region is not correctly placed, then we skip it
    // NOTE: I think we can solve this by increasing maxw
    if (prekmer.first == -1 || prekmer.second == -1 || postkmer.first == -1 ||
        postkmer.second == -1) {
      ++skip_2;
      continue;
    }
    // FIXME: understand why this is happening (chr16 on full giab genome)
    if ((uint)prekmer.second > postkmer.second + kmer_size) {
      cerr << "Error on " << qname << ". SFS starting at " << sfs.s
           << " (length " << sfs.l << ")." << endl;
    } else {
      local_extended_sfs.push_back(
          ExtSFS(string(chrom), string(qname), prekmer.second,
                 postkmer.second + kmer_size, prekmer.first,
                 postkmer.first + kmer_size));
    }
  }
  // when two SFSs are close but not overlapping, we may end up with two
  // overlapping extended SFSs merge SFS two by two until we don't need to merge
  // anything
  vector<ExtSFS> merged_extended_sfs;
  for (int i = 0; i < local_extended_sfs.size(); ++i) {
    int j;
    for (j = 0; j < merged_extended_sfs.size(); ++j) {
      if ((local_extended_sfs.at(i).s <= merged_extended_sfs.at(j).s &&
           merged_extended_sfs.at(j).s <= local_extended_sfs.at(i).e) ||
          (merged_extended_sfs.at(j).s <= local_extended_sfs.at(i).s &&
           local_extended_sfs.at(i).s <= merged_extended_sfs.at(j).e))
        break;
    }
    if (j < merged_extended_sfs.size()) {
      merged_extended_sfs[j].s =
          min(merged_extended_sfs.at(j).s, local_extended_sfs.at(i).s);
      merged_extended_sfs[j].e =
          max(merged_extended_sfs.at(j).e, local_extended_sfs.at(i).e);
      merged_extended_sfs[j].qs =
          min(merged_extended_sfs.at(j).qs, local_extended_sfs.at(i).qs);
      merged_extended_sfs[j].qe =
          max(merged_extended_sfs.at(j).qe, local_extended_sfs.at(i).qe);
    } else {
      merged_extended_sfs.push_back(local_extended_sfs.at(i));
    }
  }
  for (const auto mes : merged_extended_sfs)
    _p_extended_sfs[index].push_back(mes);

  if (lclip.second > 0) {
    _p_clips[index].push_back(
        Clip(qname, chrom, lclip.first, lclip.second, true));
  }
  if (rclip.second > 0) {
    _p_clips[index].push_back(
        Clip(qname, chrom, rclip.first, rclip.second, false));
  }
}

bool overlap(int s1, int e1, const ExtSFS &sfs) {
  int o = max(s1, sfs.s) - min(e1, sfs.e);
  return o >= 0;
}

void Extender::cluster_no_interval_tree() {
  sort(extended_sfs.begin(), extended_sfs.end());
  auto r = std::max_element(extended_sfs.begin(), extended_sfs.end(),
                            [](const ExtSFS &lhs, const ExtSFS &rhs) {
                              return lhs.e - lhs.s < rhs.e - rhs.s;
                            });
  int dist = (r->e - r->s) * 1.1;
  lprint({"Maximum extended-SFS length:", to_string(r->e - r->s),
          "bp. Using separation distance", to_string(dist) + "."});
  // find large gaps
  int prev_i = 0;
  int prev_e = extended_sfs[0].e;
  string prev_chrom = extended_sfs[0].chrom;
  vector<pair<int, int>> intervals;
  for (int i = 1; i < extended_sfs.size(); i++) {
    const auto &sfs = extended_sfs[i];
    // new chromosome
    if (sfs.chrom != prev_chrom) {
      prev_chrom = sfs.chrom;
      intervals.push_back(make_pair(prev_i, i - 1));
      prev_i = i;
      prev_e = sfs.e;
      continue;
    } else {
      if (sfs.s - prev_e > dist) {
        intervals.push_back(make_pair(prev_i, i - 1));
        prev_e = sfs.e;
        prev_i = i;
      }
    }
  }
  intervals.push_back(make_pair(prev_i, extended_sfs.size() - 1));
  // cluster each interval independently
  time_t f;
  time(&f);
  time_t s;
  time(&s);
  time_t u;
  bool printed = false;
  lprint({"Retrieved", to_string(intervals.size()), " intervals."});
  _p_sfs_clusters.resize(threads);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (int i = 0; i < intervals.size(); i++) {
    int t = omp_get_thread_num();
    interval_tree_t<int> tree;
    int j = intervals[i].first;
    int low = extended_sfs[j].s;
    int high = extended_sfs[j].e;
    int last_j = j;
    j++;
    for (; j <= intervals[i].second; j++) {
      const ExtSFS &sfs = extended_sfs[j];
      if (sfs.s <= high) {
        low = min(low, sfs.s);
        high = max(high, sfs.e);
      } else {
        for (int k = last_j; k < j;
             k++) { // CHECKME: < or <=?
                    // NOTE: <= makes the code waaaay slower
          _p_sfs_clusters[t][make_pair(low, high)].push_back(extended_sfs[k]);
        }
        low = sfs.s;
        high = sfs.e;
        last_j = j;
      }
    }
    for (int k = last_j; k <= intervals[i].second;
         k++) { // CHECKME: it was < but in that way
                // we were losing an sfs per cluster
      _p_sfs_clusters[t][make_pair(low, high)].push_back(extended_sfs[k]);
    }
    if (t == 0) {
      time(&u);
      if (u - s > 30) {
        cerr << "[I] Processed " << std::left << i
             << " intervals so far. Intervals per second: " << i / (u - f)
             << ". Time: " << u - f << "\r";
        time(&s);
        printed = true;
      }
    }
  }
  if (printed) {
    cerr << endl;
  }
}

void Extender::extract_sfs_sequences() {
  lprint({"Analyzing", to_string(ext_clusters.size()), "clusters.."});
  char *seq[threads];
  uint32_t len[threads];
  bam1_t *_p_aln[threads];
  samFile *_p_bam_file[threads];
  hts_idx_t *_p_bam_index[threads];
  bam_hdr_t *_p_bam_header[threads];
  for (int i = 0; i < threads; i++) {
    len[i] = 0;
    _p_aln[i] = bam_init1();
    _p_bam_file[i] = hts_open(config->bam.c_str(), "r");
    _p_bam_index[i] = sam_index_load(_p_bam_file[i], config->bam.c_str());
    _p_bam_header[i] = sam_hdr_read(_p_bam_file[i]);
    bgzf_mt(_p_bam_file[i]->fp.bgzf, 8, 1);
  }
  time_t f;
  time(&f);
  time_t s;
  time(&s);
  time_t u;
  bool printed = false;
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (int i = 0; i < ext_clusters.size(); i++) {
    int t = i % threads;
    const auto &cluster = ext_clusters[i];

    unordered_map<string, bool> reads;
    int cluster_s = numeric_limits<int>::max();
    int cluster_e = 0;
    for (const ExtSFS &esfs : cluster.seqs) {
      cluster_s = min(cluster_s, esfs.s);
      cluster_e = max(cluster_e, esfs.e);
      reads[esfs.qname] = true;
    }
    uint cluster_size = reads.size();
    if (cluster_size < min_w) {
      ++small_cl;
      continue;
    }

    string chrom = cluster.chrom;
    Cluster global_cluster = Cluster(chrom, cluster_s, cluster_e);

    uint cov = 0;
    string region =
        chrom + ":" + to_string(cluster_s) + "-" + to_string(cluster_e);
    hts_itr_t *itr =
        sam_itr_querys(_p_bam_index[t], _p_bam_header[t], region.c_str());
    while (sam_itr_next(_p_bam_file[t], itr, _p_aln[t]) > 0) {
      bam1_t *aln = _p_aln[t];
      if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY ||
          aln->core.flag & BAM_FSECONDARY) {
        continue;
      }
      ++cov; // FIXME: this cov takes into account also reads starting or
             // ending inside the cluster (maybe we should skip those?)

      char *qname = bam_get_qname(aln);
      if (reads.find(qname) == reads.end()) {
        continue;
      }

      uint32_t l = aln->core.l_qseq;
      if (l >= len[t]) {
        if (len[t] != 0) {
          free(seq[t]);
        }
        len[t] = l;
        seq[t] = (char *)malloc(l + 1);
      }
      uint8_t *q = bam_get_seq(aln);
      for (int i = 0; i < l; i++) {
        seq[t][i] = seq_nt16_str[bam_seqi(q, i)];
      }
      seq[t][l] = '\0';

      vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
      int qs = -1, qe = -1;
      // getting starting and ending positions on read sequence aligning to
      // start/end of cluster
      for (int i = alpairs.size() - 1; i >= 0; --i) {
        // finding starting position
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (q == -1 || r == -1) {
          continue;
        }
        if (r <= cluster_s) {
          qs = q;
          break;
        }
      }
      for (uint i = 0; i < alpairs.size(); ++i) {
        // finding ending position
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (q == -1 || r == -1) {
          continue;
        }
        if (r >= cluster_e) {
          qe = q;
          break;
        }
      }
      if (qs == -1 || qe == -1) {
        // reads starts or ends inside the cluster
        // TODO: get only remaining prefix/suffix? but this may make POA and
        // realignment harder
        ++skip_3;
      } else {
        string _seq(seq[t], qs, qe - qs + 1);
        global_cluster.add(qname, _seq);
      }
    }
    if (global_cluster.size() >= min_w) {
      global_cluster.set_cov(cov);
      _p_clusters[t].push_back(global_cluster);
      ++extcl;
    } else {
      ++small_extcl;
    }
    if (t == 0) {
      time(&u);
      if (u - s > 30) {
        cerr << "[I] Processed " << std::left << i
             << " clusters so far. Clusters per second: " << i / (u - f)
             << ". Time: " << u - f << "\r";
        time(&s);
        printed = true;
      }
    }
  }
  if (printed) {
    cerr << endl;
  }
  for (int i = 0; i < threads; i++) {
    if (len[i] > 0) {
      free(seq[i]);
    }
    bam_destroy1(_p_aln[i]);
    sam_close(_p_bam_file[i]);
  }
}

vector<Cluster> Extender::cluster_by_length(const Cluster &cluster) {
  vector<Cluster> clusters_by_len;
  for (uint c = 0; c < cluster.size(); ++c) {
    const string &name = cluster.get_name(c);
    const string &seq = cluster.get_seq(c);
    int i;
    for (i = 0; i < clusters_by_len.size(); i++) {
      float cl = clusters_by_len[i].get_len();
      float sl = seq.size();
      if (min(cl, sl) / max(cl, sl) >= config->min_ratio)
        // if (abs((int)clusters_by_len[i].get_len() - (int)seq.size()) <=
        // min_d) {
        break;
    }
    if (i == clusters_by_len.size()) {
      clusters_by_len.push_back(
          Cluster(cluster.chrom, cluster.s, cluster.e, cluster.cov));
    }
    clusters_by_len[i].add(name, seq);
  }
  return clusters_by_len;
}

vector<pair<uint, char>> Extender::parse_cigar(string cigar) {
  // -- Parsing CIGAR
  vector<pair<uint, char>> cigar_pairs;
  regex r("([0-9]+)([MIDNSHPX=])");
  regex_iterator<string::iterator> rit(cigar.begin(), cigar.end(), r);
  regex_iterator<string::iterator> rend;
  while (rit != rend) {
    int l = stoi(rit->str().substr(0, rit->str().size() - 1));
    char op = rit->str().substr(rit->str().size() - 1, 1)[0];
    cigar_pairs.push_back(make_pair(l, op));
    ++rit;
  }
  return cigar_pairs;
}

void Extender::call() {
  lprint({"Calling SVs from", to_string(clusters.size()), "clusters.."});
  time_t f;
  time(&f);
  time_t s;
  time(&s);
  time_t u;
  bool printed = false;
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (int _ = 0; _ < clusters.size(); _++) {
    int t = _ % threads;
    const Cluster &cluster = clusters[_];
    string chrom = cluster.chrom;
    if (cluster.size() < min_w) {
      continue;
    }
    const auto &clusters_by_len = cluster_by_length(cluster);
    if (config->verbose) {
      cout << "1 " << cluster.chrom << "\t" << cluster.s << "\t" << cluster.e
           << "\t" << cluster.cov << "\t" << cluster.seqs.size() << "\t";
      for (const auto s : cluster.seqs)
        cout << s.size() << ",";
      cout << "\t";
      for (const auto s : cluster.seqs)
        cout << s << ",";
      cout << endl;
      if (clusters_by_len.size() > 1) {
        for (const auto c : clusters_by_len) {
          cout << "- " << c.chrom << "\t" << c.s << "\t" << c.e << "\t" << c.cov
               << "\t" << c.seqs.size() << "\t";
          for (const auto s : c.seqs)
            cout << s.size() << ",";
          cout << "\t";
          for (const auto s : c.seqs)
            cout << s << ",";
          cout << endl;
        }
      }
    }
    // --- Sorting clusters by #sequences to get first 2 most weighted
    // clusters
    int i_max1 = -1;
    int i_max2 = -1;
    uint v_max1 = 0;
    uint v_max2 = 0;
    for (uint i = 0; i < clusters_by_len.size(); ++i) {
      if (clusters_by_len[i].size() > v_max1) {
        v_max2 = v_max1;
        i_max2 = i_max1;
        v_max1 = clusters_by_len[i].size();
        i_max1 = i;
      } else if (clusters_by_len[i].size() > v_max2) {
        v_max2 = clusters_by_len[i].size();
        i_max2 = i;
      }
    }
    vector<int> maxs({i_max1, i_max2});
    // Genotyping the two most-weighted clusters
    for (const int i : maxs) {
      if (i == -1) {
        continue;
      }
      Cluster c = clusters_by_len[i];
      if (c.size() < min_w) {
        continue;
      }
      // --- Local realignment
      vector<SV> _svs; // svs on current cluster
      string ref = string(chromosome_seqs[chrom] + c.s, c.e - c.s + 1);
      string consensus = c.poa();

      parasail_result_t *result = NULL;
      result = parasail_nw_trace_striped_16(consensus.c_str(), consensus.size(),
                                            ref.c_str(), ref.size(), 10, 1,
                                            &parasail_nuc44);
      parasail_cigar_t *cigar =
          parasail_result_get_cigar(result, consensus.c_str(), consensus.size(),
                                    ref.c_str(), ref.size(), NULL);
      string cigar_str = parasail_cigar_decode(cigar);
      if (config->verbose)
        cout << ref << "," << consensus << "," << cigar_str << "\n"
             << "#" << endl;
      int score = result->score;
      parasail_cigar_free(cigar);
      parasail_result_free(result);
      _p_alignments[t].push_back(
          Consensus(consensus, cigar_str, chrom, c.s, c.e));
      // -- Extracting SVs
      uint rpos = c.s; // position on reference
      uint cpos = 0;   // position on consensus
      auto cigar_pairs = parse_cigar(cigar_str);
      int nv = 0;
      for (const auto cigar_pair : cigar_pairs) {
        uint l = cigar_pair.first;
        char op = cigar_pair.second;
        if (op == '=' || op == 'M') {
          rpos += l;
          cpos += l;
        } else if (op == 'I') {
          if (l > config->min_sv_length) {
            SV sv = SV("INS", c.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, 1),
                       consensus.substr(cpos, l), c.size(), c.cov, nv, score,
                       false, l, cigar_str);
            sv.add_reads(c.get_names());
            _svs.push_back(sv);
            nv++;
          }
          cpos += l;
        } else if (op == 'D') {
          if (l > config->min_sv_length) {
            SV sv = SV("DEL", c.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, l),
                       string(chromosome_seqs[chrom] + rpos - 1, 1), c.size(),
                       c.cov, nv, score, false, l, cigar_str);
            sv.add_reads(c.get_names());
            _svs.push_back(sv);
            nv++;
          }
          rpos += l;
        }
      }
      for (int v = 0; v < _svs.size(); v++) {
        _svs[v].ngaps = nv;
      }
      // --- combine SVs on same consensus ---
      // vector<SV> merged_svs ;
      // std::sort(_svs.begin(), _svs.end()) ;
      // for (int i = 1; i < _svs.size(); i++) {
      //    if (_svs[i].type == "DEL" == ) {
      //        _dels.push_back(_svs[i]) ;
      //    } else {
      //        _ins.push_back(_svs[i]) ;
      //    }
      //}
      // merge only overlapping SVs
      // SV sv = _dels[0] ;
      // for (int i = 1; i < _dels.size(); i++) {
      //    auto& del = _dels[i] ;
      //    if (del.s <= sv.e) {
      //        sv.e = del.e ;
      //        int overlap = sv.e - del.s + 1 ;
      //        sv.refall = sv.altall + del..substr(overlap, del.l - overlap);
      //    }
      //}
      // -- Combine svs with same length (maybe useless now - only if diploid
      // mode)
      // vector<SV> comb_svs;
      // for (const SV &msv : merged_svs) {
      //    bool newsv_flag = true ;
      //    for (SV &csv: comb_svs) {
      //        if (abs(csv.l - msv.l) <= 10) {
      //            csv.w += msv.w ;
      //            newsv_flag = false ;
      //        }
      //    }
      //    if (newsv_flag) {
      //        comb_svs.push_back(msv);
      //    }
      //}
      for (const SV &sv : _svs) {
        _p_svs[t].push_back(sv);
      }
    }
    if (t == 0) {
      time(&u);
      if (u - s > 30) {
        cerr << "[I] Processed " << std::left << _
             << " clusters so far. Cluster per second: " << _ / (u - f)
             << ". Time: " << u - f << "\r";
        time(&s);
      }
    }
  }
  if (printed) {
    cerr << endl;
  }
}

void Extender::filter_sv_chains() {
  std::sort(svs.begin(), svs.end());
  if (config->verbose)
    for (const auto sv : svs)
      cout << sv << endl;
  if (svs.size() < 2) {
    return;
  }
  lprint({to_string(svs.size()), "SVs before chain filtering."});
  vector<SV> _svs;
  auto &prev = svs[0];
  bool reset = false;
  for (int i = 1; i < svs.size(); i++) {
    if (reset) {
      reset = false;
      prev = svs[i];
      continue;
    }
    auto &sv = svs[i];
    if (sv.chrom == prev.chrom && sv.s - prev.e < 2 * sv.l &&
        prev.type == sv.type) {
      //  check for sequence similarity
      double w_r =
          min((double)sv.w, (double)prev.w) / max((double)sv.w, (double)prev.w);
      double l_r =
          min((double)sv.l, (double)prev.l) / max((double)sv.l, (double)prev.l);
      int d = sv.s - prev.s;
      if (d < 100 && w_r >= 0.9 &&
          l_r >= config->min_ratio) { // FIXME: hardcoded + use different ratio
                                      // here. min_ratio was for clusters
        double sim;
        if (sv.type == "DEL") {
          sim = rapidfuzz::fuzz::ratio(sv.refall, prev.refall);
        } else {
          sim = rapidfuzz::fuzz::ratio(sv.altall, prev.altall);
        }
        if (sim > 70) {
          if (sv.w > prev.w) {
            _svs.push_back(sv);
          } else {
            _svs.push_back(prev);
          }
          reset = true;
          continue;
        }
      }
    }
    _svs.push_back(prev);
    prev = sv;
  }
  _svs.push_back(prev);
  svs.clear();
  svs.insert(svs.begin(), _svs.begin(), _svs.end());
  lprint({to_string(svs.size()), "SVs after chain filtering."});
}
