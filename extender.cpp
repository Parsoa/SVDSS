#include "extender.hpp"

Extender::Extender(unordered_map<string, vector<SFS>> *_SFSs) {
  SFSs = _SFSs;
  config = Configuration::getInstance();
}

void Extender::run() {
  // 0. Open BAM
  bam_file = hts_open(config->bam.c_str(), "r");
  bam_index = sam_index_load(bam_file, config->bam.c_str());
  bam_header = sam_hdr_read(bam_file);
  bgzf_mt(bam_file->fp.bgzf, 8, 1);

  // 1. Place SFSs on reference genome by extracting subalignments from read
  // alignments
  spdlog::info("Placing SFSs on reference genome");
  _p_clips.resize(config->threads);
  _p_extended_sfs.resize(config->threads);
  extend_parallel();
  for (int i = 0; i < config->threads; i++) {
    for (const auto &extsfs : _p_extended_sfs[i])
      extended_sfs.push_back(extsfs);
    clips.insert(clips.begin(), _p_clips[i].begin(), _p_clips[i].end());
  }
  spdlog::info("{}/{}/{} unplaced SFSs. {} erroneus SFSs. {} clipped SFSs.",
               unplaced, s_unplaced, e_unplaced, unknown, clips.size());

  // 2. Cluster SFSs by proximity
  spdlog::info("Clustering {} SFSs..", extended_sfs.size());
  cluster_no_interval_tree();
  map<pair<int, int>, ExtCluster> _ext_clusters;
  for (int i = 0; i < config->threads; i++)
    for (const auto &cluster : _p_sfs_clusters[i])
      _ext_clusters.insert(
          make_pair(cluster.first, ExtCluster(cluster.second)));
  for (const auto &cluster : _ext_clusters)
    ext_clusters.push_back(cluster.second);

  // 3. Extend SFSs inside each cluster to force them to start/end at same
  // reference position
  spdlog::info("Extending {} clusters..", ext_clusters.size());
  _p_clusters.resize(config->threads);
  extract_sfs_sequences();
  for (int i = 0; i < config->threads; i++)
    clusters.insert(clusters.begin(), _p_clusters[i].begin(),
                    _p_clusters[i].end());
  spdlog::info(
      "Filtered {} SFSs.Filtered {} clusters. Filtered {} global clusters.",
      unextended, small_clusters, small_extclusters);

  if (config->clusters.compare("") != 0) {
    spdlog::info("Storing clusters to {}", config->clusters);
    ofstream clofile;
    clofile.open(config->clusters);
    for (const auto &cluster : clusters) {
      // CHECKME: is ending position (cluster.first.second) inclusive? I'm
      // assuming yes
      clofile << cluster.chrom << ":" << cluster.s + 1 << "-" << cluster.e + 1
              << "\t" << cluster.size();
      for (size_t i = 0; i < cluster.size(); ++i)
        clofile << "\t" << cluster.get_name(i) << ":" << cluster.get_seq(i);
      clofile << endl;
    }
    clofile.close();
  }

  // 4. Call SVs
  spdlog::info("Calling SVs from {} clusters..", clusters.size());
  _p_svs.resize(config->threads);
  _p_alignments.resize(config->threads);
  call();
  for (int i = 0; i < config->threads; i++) {
    svs.insert(svs.begin(), _p_svs[i].begin(), _p_svs[i].end());
    alignments.insert(alignments.begin(), _p_alignments[i].begin(),
                      _p_alignments[i].end());
  }
  filter_sv_chains();
  sam_close(bam_file);
  spdlog::info("Total SVs: {}", svs.size());
}

/* Get first/last kmer entirely mapped and with a single occurrence in a
 * subalignment of interest. Return the kmer as the initial pair of positions in
 * the alignment */
pair<int, int> Extender::get_unique_kmers(const vector<pair<int, int>> &alpairs,
                                          const uint k, const bool from_end,
                                          string chrom) {
  if (alpairs.size() < k)
    return make_pair(-1, -1);

  // Do kmer counting in region
  map<string, int> kmers;
  string kmer_seq;
  size_t i = 0;
  while (i < alpairs.size() - k + 1) {
    bool skip = false;
    for (size_t j = i; j < i + k; j++) {
      if (alpairs[j].first == -1 || alpairs[j].second == -1) {
        // we want clean kmers only - ie placed kmers, no insertions or
        // deletions
        skip = true;
        i = j + 1; // jump to next clean/placed position
        break;
      }
    }
    if (skip)
      continue;
    string kmer(chromosome_seqs[chrom] + alpairs[i].second, k);
    ++kmers[kmer];
    ++i;
  }
  // Get first/last kmer with single occurrence
  pair<int, int> last_kmer = make_pair(-1, -1);
  i = 0;
  while (i < alpairs.size() - k + 1) {
    int offset = i;
    if (from_end)
      offset = alpairs.size() - k - i;
    assert(offset >= 0);
    bool skip = false;
    for (size_t j = offset; j < offset + k; j++) {
      if (alpairs[j].first == -1 || alpairs[j].second == -1) {
        skip = true;
        i += (j - offset); // jump to next possible start position
        break;
      }
    }
    if (skip) {
      ++i;
      continue;
    }
    last_kmer = alpairs[offset];
    string kmer(chromosome_seqs[chrom] + alpairs[offset].second, k);
    if (kmers[kmer] == 1)
      break;
    ++i;
  }
  return last_kmer;
}

// Parallelize within each chromosome
void Extender::extend_parallel() {
  // Initializing
  for (int i = 0; i < 2; i++) {
    bam_entries.push_back(vector<vector<bam1_t *>>(config->threads));
    for (int j = 0; j < config->threads; j++)
      for (int k = 0; k < config->batch_size / config->threads; k++)
        bam_entries[i][j].push_back(bam_init1());
  }

  int p = 0;
  int b = 0;
  load_batch_bam(p);
  time_t start_time;
  time_t curr_time;
  time(&start_time);
  bool should_load = true;
  bool should_process = true;
  while (should_process) {
    if (!should_load)
      should_process = false;
#pragma omp parallel for num_threads(config->threads + 1)
    // TODO: avoid additional thread. Make -1 before (same in other classes,
    // with -2)
    for (int i = 0; i < config->threads + 1; i++) {
      int t = omp_get_thread_num();
      if (t == 0) {
        // First thread loads next batch
        if (should_load)
          should_load = load_batch_bam((p + 1) % 2);
      } else
        // Other threads process batch
        extend_batch(bam_entries[p][t - 1], t - 1);
    }
    p += 1;
    p %= 2;
    b += 1;
    time(&curr_time);
    if (curr_time - start_time == 0)
      ++curr_time;
    cerr << "Extended batch " << b << ". Time: " << curr_time - start_time
         << "\r";
  }

  // cleanup
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < config->threads; j++)
      for (int k = 0; k < config->batch_size / config->threads; k++)
        bam_destroy1(bam_entries[i][j][k]);
}

/* Load batch from BAM file and store to input entry p. The logic behind is:
 * fill position i per each thread, then move to position i+1.. */
bool Extender::load_batch_bam(int p) {
  int i = 0;
  int nseqs = 0;
  while (sam_read1(bam_file, bam_header,
                   bam_entries[p][nseqs % config->threads][i]) >= 0) {
    bam1_t *aln = bam_entries[p][nseqs % config->threads][i];
    if (aln == nullptr) {
      spdlog::critical("nullptr. Why are we here? Please check");
      exit(1);
    }
    if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY ||
        aln->core.flag & BAM_FSECONDARY) {
      spdlog::warn("Non primary alignment.. Is bam file smoothed?");
      continue;
    }
    char *qname = bam_get_qname(aln);
    if (SFSs->find(qname) == SFSs->end())
      continue;
    ++nseqs;
    if (nseqs % config->threads == 0)
      ++i;
    if (nseqs == config->batch_size)
      return true;
  }
  // last batch is incomplete since we reached the end of .bam file
  // TODO: can we do like in ping_pong?
  if (nseqs != config->batch_size) {
    for (int j = nseqs % config->threads; j < config->threads; j++)
      for (int _ = i; _ < config->batch_size / config->threads; _++)
        bam_entries[p][j][_] = nullptr;
    for (int j = 0; j < nseqs % config->threads; j++)
      for (int _ = i + 1; _ < config->batch_size / config->threads; _++)
        bam_entries[p][j][_] = nullptr;
  }
  return false;
}

/* Extend a batch of alignments */
void Extender::extend_batch(vector<bam1_t *> bam_entries, int index) {
  bam1_t *aln;
  for (size_t b = 0; b < bam_entries.size(); b++) {
    aln = bam_entries[b];
    if (aln == nullptr)
      break;
    extend_alignment(aln, index);
  }
}

/* Place SFSs on the alignment and extend them using k-mers */
void Extender::extend_alignment(bam1_t *aln, int index) {
  char *qname = bam_get_qname(aln);
  uint32_t *cigar = bam_get_cigar(aln);
  vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
  string chrom(bam_header->target_name[aln->core.tid]);
  // TODO: do this check while loading the alignments
  if (chromosome_seqs.find(chrom) == chromosome_seqs.end())
    return;

  // NOTE: we may have more sfs on a clipped read, but all of them will produce
  // the same clip
  pair<uint, uint> lclip = make_pair(0, 0);
  pair<uint, uint> rclip = make_pair(0, 0);
  int last_pos = 0;
  vector<ExtSFS> local_extended_sfs;
  for (const SFS &sfs : SFSs->at(qname)) {
    int s = sfs.s;
    int e = sfs.s + sfs.l - 1;
    int hp_tag = sfs.htag;
    int aln_start = -1;
    int aln_end = -1;
    vector<pair<int, int>> local_alpairs; // subalignment
    // find start and end of SFS in read's alignment
    int refs = -1;
    int refe = -1;
    for (size_t i = last_pos; i < alpairs.size(); i++) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      if (q == -1 || r == -1)
        continue;
      else if (q < s) {
        // here we are getting the last "placed" base before the SFS
        // <= seems more correct to me but using < we are more flexible
        last_pos = i;
        refs = r;
        aln_start = i;
      } else if (q > e) {
        // here we are getting the first "placed" base after the SFS
        // >= seems more correct but using > we are more flexible
        refe = r;
        aln_end = i;
        break;
      }
    }
    // Current SFS is aligned from refs to refe. In case of an insertion, the
    // interval will cover the entire insertion

    // We extract the local alignment of the region of interest
    if (refs == -1 && refe == -1) {
      // we couldn't place the first and the last base, so we skip this -
      // otherwise we'll end up considering the entire read
      ++unplaced;
      continue;
    } else if (refs == -1) {
      uint op = bam_cigar_op(*(cigar + 0));
      uint l = bam_cigar_oplen(*(cigar + 0));
      if (op == BAM_CSOFT_CLIP && config->clipped)
        lclip = make_pair(aln->core.pos, l);
      else
        ++s_unplaced;
      continue; // in any case, we skip this SFS
    } else if (refe == -1) {
      uint op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
      uint l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
      if (op == BAM_CSOFT_CLIP && config->clipped)
        rclip = make_pair(bam_endpos(aln), l);
      else
        ++e_unplaced;
      continue; // in any case, we skip this SFS
    } else {
      // we placed the first and last base, so we extract the subalignment
      int last_r = refs - 1;
      for (int i = aln_start; i <= aln_end; i++) {
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (r == -1) {
          if (refs <= last_r && last_r <= refe)
            local_alpairs.push_back(make_pair(q, r));
        } else {
          last_r = r;
          if (refs <= r && r <= refe)
            local_alpairs.push_back(make_pair(q, r));
        }
        // We break when we found a placed base at or after the reference end
        if (q != -1 && r != -1 && r >= refe)
          break;
      }
    }

    // SFS has been placed and local_alpairs contains the subalignment

    // extract the config->flank pairs preceding the region of interest
    vector<pair<int, int>> pre_alpairs;
    uint n = 0;
    for (int i = aln_start - 1; i >= 0; --i) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      pre_alpairs.push_back(make_pair(q, r));
      ++n;
      if (n == config->flank)
        break;
    }
    reverse(pre_alpairs.begin(), pre_alpairs.end());
    // extract the config->flank pairs following the region of interest
    vector<pair<int, int>> post_alpairs;
    n = 0;
    for (uint i = aln_end + 1; i < alpairs.size(); i++) {
      int q = alpairs[i].first;
      int r = alpairs[i].second;
      post_alpairs.push_back(make_pair(q, r));
      ++n;
      if (n == config->flank)
        break;
    }

    // get the unique kmer in the upstream and downstream config->flank bp
    // regions
    pair<int, int> prekmer =
        get_unique_kmers(pre_alpairs, config->ksize, true,
                         chrom); // true for first kmer found (shorter cluster)
    pair<int, int> postkmer =
        get_unique_kmers(post_alpairs, config->ksize, false,
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
    // NOTE: I think we can solve this by increasing config->flank.
    if (prekmer.first == -1 || prekmer.second == -1 || postkmer.first == -1 ||
        postkmer.second == -1) {
      spdlog::warn("SFS has not been placed. But why? Check this plz.");
      ++unknown;
      continue;
    }
    // FIXME: understand why this is happening (chr16 on full giab genome)
    if ((uint)prekmer.second > postkmer.second + config->ksize) {
      spdlog::warn("Error on {}. SFS starting at {} (length {})", qname, sfs.s,
                   sfs.l);
    } else {
      local_extended_sfs.push_back(
          ExtSFS(string(chrom), string(qname), prekmer.second,
                 postkmer.second + config->ksize, prekmer.first,
                 postkmer.first + config->ksize, hp_tag));
    }
  }
  // When two SFSs are close but not overlapping, we may end up with two
  // overlapping extended SFSs. We need to merge SFSs two by two until we don't
  // need to merge anything
  vector<ExtSFS> merged_extended_sfs;
  for (size_t i = 0; i < local_extended_sfs.size(); ++i) {
    size_t j;
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

  // if (lclip.second > 0)
  //   _p_clips[index].push_back(
  //       Clip(qname, chrom, lclip.first, lclip.second, true));
  // if (rclip.second > 0)
  //   _p_clips[index].push_back(
  //       Clip(qname, chrom, rclip.first, rclip.second, false));
}

void Extender::cluster_no_interval_tree() {
  sort(extended_sfs.begin(), extended_sfs.end());
  auto r = max_element(extended_sfs.begin(), extended_sfs.end(),
                       [](const ExtSFS &lhs, const ExtSFS &rhs) {
                         return lhs.e - lhs.s < rhs.e - rhs.s;
                       });
  int dist = (r->e - r->s) * 1.1;
  spdlog::info(
      "Maximum extended SFS length: {}bp. Using separation distance {}.",
      r->e - r->s, dist);
  // find large gaps
  int prev_i = 0;
  int prev_e = extended_sfs[0].e;
  string prev_chrom = extended_sfs[0].chrom;
  vector<pair<int, int>> intervals;
  for (size_t i = 1; i < extended_sfs.size(); i++) {
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
  _p_sfs_clusters.resize(config->threads);
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < intervals.size(); i++) {
    int t = omp_get_thread_num();
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
  }
}

/* Assign coverage and read (sub)sequence to each cluster  */
void Extender::extract_sfs_sequences() {
  // Allocate
  char *seq[config->threads];
  uint32_t len[config->threads];
  bam1_t *_p_aln[config->threads];
  samFile *_p_bam_file[config->threads];
  hts_idx_t *_p_bam_index[config->threads];
  bam_hdr_t *_p_bam_header[config->threads];
  for (int i = 0; i < config->threads; i++) {
    len[i] = 0;
    _p_aln[i] = bam_init1();
    _p_bam_file[i] = hts_open(config->bam.c_str(), "r");
    _p_bam_index[i] = sam_index_load(_p_bam_file[i], config->bam.c_str());
    _p_bam_header[i] = sam_hdr_read(_p_bam_file[i]);
    bgzf_mt(_p_bam_file[i]->fp.bgzf, 8, 1);
  }

#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < ext_clusters.size(); i++) {
    int t = omp_get_thread_num();
    const auto &cluster = ext_clusters[i];

    // Force all extended SFSs to start and end at the same position. Build a
    // "global" cluster
    set<string> reads;
    int cluster_s = numeric_limits<int>::max();
    int cluster_e = 0;
    for (const ExtSFS &esfs : cluster.seqs) {
      cluster_s = min(cluster_s, esfs.s);
      cluster_e = max(cluster_e, esfs.e);
      reads.insert(esfs.qname);
    }

    size_t cluster_size = reads.size();
    if (cluster_size < config->min_cluster_weight) {
      ++small_clusters;
      continue;
    }

    Cluster global_cluster = Cluster(cluster.chrom, cluster_s, cluster_e);

    // Iterate over alignments falling in the cluster region to: (i) get total
    // number of reads and (ii) get SFS sequence, one per read
    uint cov = 0;
    string region =
        cluster.chrom + ":" + to_string(cluster_s) + "-" + to_string(cluster_e);
    hts_itr_t *itr =
        sam_itr_querys(_p_bam_index[t], _p_bam_header[t], region.c_str());
    while (sam_itr_next(_p_bam_file[t], itr, _p_aln[t]) > 0) {
      bam1_t *aln = _p_aln[t];
      if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY ||
          aln->core.flag & BAM_FSECONDARY)
        continue;
      ++cov; // FIXME: this cov takes into account also reads starting or
             // ending inside the cluster (maybe we should skip those?)

      char *qname = bam_get_qname(aln);
      if (reads.find(qname) == reads.end())
        continue;

      // If we have a SFS on this read, load the read sequence
      uint32_t l = aln->core.l_qseq;
      if (l >= len[t]) {
        if (len[t] != 0)
          free(seq[t]);
        len[t] = l;
        seq[t] = (char *)malloc(l + 1);
      }
      uint8_t *q = bam_get_seq(aln);
      for (uint i = 0; i < l; i++)
        seq[t][i] = seq_nt16_str[bam_seqi(q, i)];
      seq[t][l] = '\0';

      // Extract "global" sequence by getting start/end position on read
      // sequence aligning to start/end of cluster
      vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
      int qs = -1, qe = -1;
      for (int i = alpairs.size() - 1; i >= 0; --i) {
        // finding starting position
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (q == -1 || r == -1)
          continue;
        if (r <= cluster_s) {
          qs = q;
          break;
        }
      }
      for (uint i = 0; i < alpairs.size(); ++i) {
        // finding ending position
        int q = alpairs[i].first;
        int r = alpairs[i].second;
        if (q == -1 || r == -1)
          continue;
        if (r >= cluster_e) {
          qe = q;
          break;
        }
      }
      if (qs == -1 || qe == -1) {
        // reads starts or ends inside the cluster
        // TODO: get only remaining prefix/suffix? but this may make POA and
        // realignment harder
        ++unextended;
      } else {
        string _seq(seq[t], qs, qe - qs + 1);
        global_cluster.add(qname, _seq);
      }
    }
    if (global_cluster.size() >= config->min_cluster_weight) {
      global_cluster.set_cov(cov);
      _p_clusters[t].push_back(global_cluster);
    } else
      ++small_extclusters;
  }

  // clean
  for (int i = 0; i < config->threads; i++) {
    if (len[i] > 0)
      free(seq[i]);
    bam_destroy1(_p_aln[i]);
    sam_close(_p_bam_file[i]);
  }
}

/* Split cluster in subclusters */
vector<Cluster> Extender::cluster_by_length(const Cluster &cluster) {
  vector<Cluster> clusters_by_len;
  for (uint c = 0; c < cluster.size(); ++c) {
    const string &name = cluster.get_name(c);
    const string &seq = cluster.get_seq(c);
    size_t i;
    for (i = 0; i < clusters_by_len.size(); i++) {
      float cl = clusters_by_len[i].get_len();
      float sl = seq.size();
      if (min(cl, sl) / max(cl, sl) >= config->min_ratio)
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

/* Convert a CIGAR string into a vector of pairs */
vector<pair<uint, char>> Extender::parse_cigar(string cigar) {
  // TODO: we already have a parse_cigar in bam.hpp
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

/* Call SVs by POA+realignment */
void Extender::call() {
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < clusters.size(); i++) {
    int t = omp_get_thread_num();
    const Cluster &cluster = clusters[i];
    string chrom = cluster.chrom;
    const auto &clusters_by_len = cluster_by_length(cluster);

    // Sorting clusters by #sequences to get first 2 most weighted clusters
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
      if (c.size() < config->min_cluster_weight)
        continue;

      vector<SV> _svs;

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
                       string(chromosome_seqs[chrom] + rpos - 1, 1) +
                           consensus.substr(cpos, l),
                       c.size(), c.cov, nv, score, false, l, cigar_str);
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
      for (size_t v = 0; v < _svs.size(); v++)
        _svs[v].ngaps = nv;
      for (const SV &sv : _svs)
        _p_svs[t].push_back(sv);
    }
  }
}

/* Merge close and similar SVs */
void Extender::filter_sv_chains() {
  if (svs.size() < 2)
    return;
  sort(svs.begin(), svs.end());
  spdlog::info("{} SVs before chain filtering.", svs.size());
  vector<SV> _svs;
  auto &prev = svs[0];
  bool reset = false;
  for (size_t i = 1; i < svs.size(); i++) {
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
        if (sv.type == "DEL")
          sim = rapidfuzz::fuzz::ratio(sv.refall, prev.refall);
        else
          sim = rapidfuzz::fuzz::ratio(sv.altall, prev.altall);
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
}
