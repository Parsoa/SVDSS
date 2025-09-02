#include "caller.hpp"

void Caller::run() {
  config = Configuration::getInstance();

  // load reference genome and SFS
  load_chromosomes(config->reference);

  SFSs = parse_sfsfile(config->sfs);

  Clusterer C = Clusterer(&SFSs);
  C.run();

  spdlog::info("Calling SVs from {} clusters..", C.clusters.size());
  _p_svs.resize(config->threads);
  _p_alignments.resize(config->threads);
  pcall(C.clusters);
  for (int i = 0; i < config->threads; i++) {
    svs.insert(svs.begin(), _p_svs[i].begin(), _p_svs[i].end());
    alignments.insert(alignments.begin(), _p_alignments[i].begin(),
                      _p_alignments[i].end());
  }
  sort(svs.begin(), svs.end());
  clean_dups();
  spdlog::info("{} SVs before chain filtering.", svs.size());
  filter_sv_chains();
  spdlog::info("Writing {} SVs.", svs.size());
  sort(svs.begin(), svs.end());
  write_vcf();

  if (config->poa.compare("") != 0) {
    spdlog::info("Writing POA alignments to {}..", config->poa);
    write_sam();
  }

  if (config->clipped) {
    spdlog::warn(
        "Calling imprecise SVs from clipped alignments is experimental");
    interval_tree_t<int> vartree;
    for (const auto &sv : svs)
      vartree.insert({sv.s - 1000, sv.e + 1000});
    vector<SV> clipped_svs;
    Clipper clipper(C.clips);
    clipper.call(config->threads, vartree);
    int s = 0;
    for (int i = 0; i < config->threads; i++) {
      s += clipper._p_svs[i].size();
      clipped_svs.insert(clipped_svs.begin(), clipper._p_svs[i].begin(),
                         clipper._p_svs[i].end());
    }
    spdlog::info("Predicted {} SVs from clipped alignments", s);
    for (const SV &sv : clipped_svs)
      cout << sv << endl;
  }

  destroy_chromosomes();
}

void Caller::write_vcf() {
  print_vcf_header();
  for (const SV &sv : svs)
    cout << sv << endl;
}

void Caller::write_sam() {
  ofstream osam;
  osam.open(config->poa);
  osam << "@HD\tVN:1.4" << endl;
  for (size_t i = 0; i < chromosomes.size(); ++i)
    osam << "@SQ\tSN:" << chromosomes[i] << "\t"
         << "LN:" << strlen(chromosome_seqs[chromosomes[i]]) << endl;
  for (size_t j = 0; j < alignments.size(); j++)
    osam << alignments[j] << endl;
  osam.close();
}

// /* Split cluster in subclusters based on length */
vector<Cluster> Caller::split_cluster_by_len(const Cluster &cluster) {
  vector<Cluster> subclusters;
  for (uint c = 0; c < cluster.size(); ++c) {
    const SubRead &sr = cluster.get_subread(c);
    size_t i;
    for (i = 0; i < subclusters.size(); i++) {
      float cl = subclusters[i].get_len();
      float sl = sr.size();
      if (min(cl, sl) / max(cl, sl) >= config->min_ratio)
        break;
    }
    if (i == subclusters.size()) {
      subclusters.push_back(Cluster(cluster.chrom, cluster.s, cluster.e,
                                    cluster.cov, cluster.cov0, cluster.cov1,
                                    cluster.cov2));
    }
    subclusters[i].add_subread(sr);
  }
  return subclusters;
}

// Split cluster in subclusters
vector<Cluster> Caller::split_cluster(const Cluster &cluster) {
  // Step 1: split cluster by haplotype tag
  Cluster cluster_0 = cluster;
  cluster_0.clear();
  Cluster cluster_1 = cluster;
  cluster_1.clear();
  Cluster cluster_2 = cluster;
  cluster_2.clear();
  for (const SubRead &sr : cluster.subreads) {
    if (config->useht) {
      if (sr.htag == 1)
        cluster_1.add_subread(sr);
      else if (sr.htag == 2)
        cluster_2.add_subread(sr);
      else
        // 0 or no tag
        cluster_0.add_subread(sr);
    } else {
      cluster_0.add_subread(sr);
    }
  }
  cluster_0.cov1 = -1;
  cluster_0.cov2 = -1;
  cluster_1.cov0 = -1;
  cluster_1.cov2 = -1;
  cluster_2.cov0 = -1;
  cluster_2.cov1 = -1;

  vector<Cluster> out_subclusters;
  if (cluster_1.size() == 0 && cluster_2.size() == 0) {
    // no alignment is tagged, use length
    vector<Cluster> subclusters = split_cluster_by_len(cluster_0);
    int i_max1 = -1, i_max2 = -1;
    uint v_max1 = 0, v_max2 = 0;
    for (uint i = 0; i < subclusters.size(); ++i) {
      if (subclusters[i].size() > v_max1) {
        v_max2 = v_max1;
        i_max2 = i_max1;
        v_max1 = subclusters[i].size();
        i_max1 = i;
      } else if (subclusters[i].size() > v_max2) {
        v_max2 = subclusters[i].size();
        i_max2 = i;
      }
    }
    if (i_max1 != -1)
      out_subclusters.push_back(subclusters[i_max1]);
    if (i_max2 != -1)
      out_subclusters.push_back(subclusters[i_max2]);
  } else {
    int both = (cluster_1.size() > 0 ? 1 : 0) + (cluster_2.size() > 0 ? 2 : 0);
    vector<Cluster> subclusters_1 = split_cluster_by_len(cluster_1);
    vector<Cluster> subclusters_2 = split_cluster_by_len(cluster_2);
    Cluster new_cluster(cluster.chrom, cluster.s, cluster.e, cluster.cov,
                        cluster.cov0, -1, -1);

    for (uint c = 0; c < cluster_0.size(); ++c) {
      const SubRead &sr = cluster_0.get_subread(c);
      float sl = sr.size();

      // check if we have to put this subread in 1 or 2
      int best_1 = -1;
      int best_ratio_1 = -1;
      for (uint i = 0; i < subclusters_1.size(); i++) {
        float cl = subclusters_1[i].get_len();
        float r = min(cl, sl) / max(cl, sl);
        if (r >= config->min_ratio && r > best_ratio_1) {
          best_1 = i;
          best_ratio_1 = r;
        }
      }
      int best_2 = -1;
      int best_ratio_2 = -1;
      for (uint i = 0; i < subclusters_2.size(); i++) {
        float cl = subclusters_2[i].get_len();
        float r = min(cl, sl) / max(cl, sl);
        if (r >= config->min_ratio && r > best_ratio_2) {
          best_2 = i;
          best_ratio_2 = r;
        }
      }

      if (both == 1) {
        assert(best_2 == -1);
        if (best_1 == -1)
          new_cluster.add_subread(sr);
        else {
          subclusters_1[best_1].add_subread(sr);
          ++subclusters_1[best_1].cov1;
          --new_cluster.cov0;
        }
      } else if (both == 2) {
        assert(best_1 == -1);
        if (best_2 == -1)
          new_cluster.add_subread(sr);
        else {
          subclusters_2[best_2].add_subread(sr);
          ++subclusters_2[best_2].cov2;
          --new_cluster.cov0;
        }
      } else {
        if (best_1 != -1 && best_ratio_1 > best_ratio_2) {
          subclusters_1[best_1].add_subread(sr);
          ++subclusters_1[best_1].cov1;
          --new_cluster.cov0;
        } else if (best_2 != -1 && best_ratio_2 > best_ratio_1) {
          subclusters_2[best_2].add_subread(sr);
          ++subclusters_2[best_2].cov2;
          --new_cluster.cov0;
        } else {
        }
      }
    }

    uint v_max = 0;
    int i_max = -1;
    for (uint i = 0; i < subclusters_1.size(); ++i) {
      if (subclusters_1[i].size() > v_max) {
        v_max = subclusters_1[i].size();
        i_max = i;
      }
    }
    if (i_max != -1)
      out_subclusters.push_back(subclusters_1[i_max]);
    v_max = 0, i_max = -1;
    for (uint i = 0; i < subclusters_2.size(); ++i) {
      if (subclusters_2[i].size() > v_max) {
        v_max = subclusters_2[i].size();
        i_max = i;
      }
    }
    if (i_max != -1)
      out_subclusters.push_back(subclusters_2[i_max]);

    if (both != 3) {
      vector<Cluster> new_subclusters = split_cluster_by_len(new_cluster);
      v_max = 0, i_max = -1;
      for (uint i = 0; i < new_subclusters.size(); ++i) {
        if (new_subclusters[i].size() > v_max) {
          v_max = new_subclusters[i].size();
          i_max = i;
        }
      }
      if (i_max != -1) {
        if (both == 1)
          new_subclusters[i_max].cov1 = -1;
        else
          new_subclusters[i_max].cov2 = -1;
        out_subclusters.push_back(new_subclusters[i_max]);
      }
    }
  }

  assert(out_subclusters.size() > 0 && out_subclusters.size() <= 2);
  return out_subclusters;
}

string Caller::run_poa(const vector<string> &seqs) {
  uint n_seqs = seqs.size();
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->align_mode = 0; // global
  abpt->disable_seeding = 1;
  abpt->progressive_poa = 0;
  abpt->amb_strand = 0;
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  // abpt->is_diploid = 1; // TODO: maybe this works now
  // abpt->max_n_cons = 2;
  // abpt->min_freq = 0.25;
  abpoa_post_set_para(abpt);

  // abpt->match = 2;      // match score
  // abpt->mismatch = 4;   // mismatch penalty
  // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
  // abpt->gap_open1 = 4;  // gap open penalty #1
  // abpt->gap_ext1 = 2;   // gap extension penalty #1
  // abpt->gap_open2 = 24; // gap open penalty #2
  // abpt->gap_ext2 = 1;   // gap extension penalty #2
  // gap_penalty = min{gap_open1 + gap_len*gap_ext1, gap_open2+gap_len*gap_ext2}

  int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
  uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
  for (uint i = 0; i < n_seqs; ++i) {
    seq_lens[i] = seqs[i].size();
    bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
    for (int j = 0; j < seq_lens[i]; ++j)
      bseqs[i][j] = _char26_table[(int)seqs[i][j]];
  }

  abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);
  abpoa_cons_t *abc = ab->abc;
  string cons = ""; // XXX: we may avoid converting to ACGT here since we need
                    // to reconvert back for ksw2
  if (abc->n_cons > 0)
    for (int j = 0; j < abc->cons_len[0]; ++j)
      cons += "ACGTN"[abc->cons_base[0][j]];

  for (uint i = 0; i < n_seqs; ++i)
    free(bseqs[i]);

  free(bseqs);
  free(seq_lens);
  abpoa_free(ab);
  abpoa_free_para(abpt);

  return cons;
}

// Call SVs by POA+realignment
void Caller::pcall(const vector<Cluster> &clusters) {
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < clusters.size(); i++) {
    int t = omp_get_thread_num();
    const Cluster &cluster = clusters[i];
    if (cluster.size() < config->min_cluster_weight)
      continue;
    string chrom = cluster.chrom;

    const vector<Cluster> &subclusters = split_cluster(cluster);

    // Calling from one or two clusters
    for (const Cluster &cl : subclusters) {
      // if (cl.size() < config->min_cluster_weight)
      //   continue;

      vector<SV> _svs;

      string ref = string(chromosome_seqs[chrom] + cl.s, cl.e - cl.s + 1);
      string consensus = run_poa(cl.get_seqs());

      // ksw2 stuff - TODO: move to a separate function
      int sc_mch = 1, sc_mis = -9, gapo = 16, gape = 2, gapo2 = 41, gape2 = 1;
      int8_t a = (int8_t)sc_mch,
             b = sc_mis < 0 ? (int8_t)sc_mis : -(int8_t)sc_mis; // a>0 and b<0
      int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                        b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
      uint tl = ref.size(), ql = consensus.size();
      uint8_t *ts = (uint8_t *)malloc(tl);
      uint8_t *qs = (uint8_t *)malloc(ql);
      for (i = 0; i < tl; ++i)
        ts[i] = _char26_table[(uint8_t)ref[i]]; // encode to 0/1/2/3
      for (i = 0; i < ql; ++i)
        qs[i] = _char26_table[(uint8_t)consensus[i]];

      ksw_extz_t ez;
      memset(&ez, 0, sizeof(ksw_extz_t));
      ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo2, gape2, -1, -1,
                    -1, 0, &ez);

      int score = ez.score;
      string cigar_str = "";
      for (int i = 0; i < ez.n_cigar; ++i) {
        cigar_str += to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];
      }
      _p_alignments[t].push_back(
          Consensus(consensus, cigar_str, chrom, cl.s, cl.e));

      // -- Extracting SVs
      uint rpos = cl.s; // position on reference
      uint cpos = 0;    // position on consensus
      CIGAR cigar;
      cigar.parse_cigar(cigar_str.c_str());
      int nv = 0;
      for (const auto &cigar_pair : cigar.ops) {
        uint l = cigar_pair.first;
        char op = cigar_pair.second;
        if (op == '=' || op == 'M') {
          rpos += l;
          cpos += l;
        } else if (op == 'I') {
          if (l >= config->min_sv_length) {
            SV sv = SV("INS", cl.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, 1),
                       string(chromosome_seqs[chrom] + rpos - 1, 1) +
                           consensus.substr(cpos, l),
                       cl.size(), cl.cov, nv, score, false, l, cigar_str);
            sv.add_reads(cl.get_names());
            _svs.push_back(sv);
            nv++;
          }
          cpos += l;
        } else if (op == 'D') {
          if (l >= config->min_sv_length) {
            SV sv = SV("DEL", cl.chrom, rpos,
                       string(chromosome_seqs[chrom] + rpos - 1, l),
                       string(chromosome_seqs[chrom] + rpos - 1, 1), cl.size(),
                       cl.cov, nv, score, false, l, cigar_str);
            sv.add_reads(cl.get_names());
            _svs.push_back(sv);
            nv++;
          }
          rpos += l;
        }
      }
      for (size_t v = 0; v < _svs.size(); v++) {
        _svs[v].ngaps = nv;
        _svs[v].set_gt("./.", 100);
        _svs[v].set_cov(cl.cov, cl.cov0, cl.cov1, cl.cov2);
        _svs[v].set_rvec(cluster.reads);
      }
      for (const SV &sv : _svs)
        _p_svs[t].push_back(sv);
    }
  }
}

// Clean same SV reported twice
void Caller::clean_dups() {
  vector<SV> _svs;
  string last_chrom = "";
  int last_pos = -1;
  string last_refall = "";
  string last_altall = "";
  for (size_t i = 0; i < svs.size(); i++) {
    if (last_chrom != svs[i].chrom || last_pos != svs[i].s ||
        last_refall != svs[i].refall || last_altall != svs[i].altall)
      _svs.push_back(svs[i]);
    last_chrom = svs[i].chrom;
    last_pos = svs[i].s;
    last_refall = svs[i].refall;
    last_altall = svs[i].altall;
  }
  svs.clear();
  svs.insert(svs.begin(), _svs.begin(), _svs.end());
}

/* Merge close and similar SVs */
void Caller::filter_sv_chains() {
  if (svs.size() < 2)
    return;

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
        if (sim > 70) { // FIXME: hardcoded
          if (sv.w > prev.w)
            _svs.push_back(sv);
          else
            _svs.push_back(prev);
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

void Caller::print_vcf_header() {
  cout << "##fileformat=VCFv4.2" << endl;
  cout << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
          "data_collections/HGSVC2/technical/reference/"
          "20200513_hg38_NoALT/"
          "hg38.no_alt.fa.gz"
       << endl;
  for (size_t i = 0; i < chromosomes.size(); ++i) {
    cout << "##contig=<ID=" << chromosomes[i]
         << ",length=" << strlen(chromosome_seqs[chromosomes[i]]) << ">"
         << endl;
  }
  cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
  cout << "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant "
          "class\">"
       << endl;
  cout << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant "
          "type\">"
       << endl;
  cout << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="
          "\"Difference in "
          "length between REF and ALT alleles\">"
       << endl;
  cout << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End "
          "position of "
          "the variant described in this record\">"
       << endl;
  cout << "##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number "
          "of "
          "alignments supporting this record\">"
       << endl;
  cout << "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus\">"
       << endl;
  cout << "##INFO=<ID=COV0,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus (no HP)\">"
       << endl;
  cout << "##INFO=<ID=COV1,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus (HP=1)\">"
       << endl;
  cout << "##INFO=<ID=COV2,Number=1,Type=Integer,Description=\"Total "
          "number of "
          "alignments covering this locus (HP=2)\">"
       << endl;
  cout << "##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Alignment "
          "score\">"
       << endl;
  cout << "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of "
          "variations on same consensus\">"
       << endl;
  cout << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="
          "\"Imprecise "
          "structural variation\">"
       << endl;
  cout << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR of "
          "consensus\">"
       << endl;
  cout << "##INFO=<ID=READS,Number=.,Type=String,Description=\"Reads "
          "identifiers supporting the call\">"
       << endl;
  cout << "##INFO=<ID=RVEC,Number=.,Type=String,Description=\"Reads vector "
          "used by genotyper\">"
       << endl;
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
       << endl;
  cout << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype "
          "quality\">"
       << endl;
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEFAULT"
       << endl;
}
