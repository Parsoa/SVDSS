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
    spdlog::critical(
        "Imprecise SVs from clipped alignments not currently supported");
    // interval_tree_t<int> vartree;
    // for (const auto &sv : svs)
    //   vartree.insert({sv.s - 1000, sv.e + 1000});
    // vector<SV> clipped_svs;
    // Clipper clipper(extender.clips);
    // clipper.call(config->threads, vartree);
    // int s = 0;
    // for (int i = 0; i < config->threads; i++) {
    //   s += clipper._p_svs[i].size();
    //   clipped_svs.insert(svs.begin(), clipper._p_svs[i].begin(),
    //                      clipper._p_svs[i].end());
    // }
    // cerr << "Predicted " << s << " SVs from clipped SFS." << endl;
    // string vcf_path = config->workdir + "/svs_clipped.vcf";
    // cerr << "Exporting " << clipped_svs.size() << " SV calls to " << vcf_path
    //      << ".." << endl;
    // for (const SV &sv : clipped_svs)
    //   cout << sv << endl;
  }
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

/* Split cluster in subclusters */
vector<Cluster> Caller::split_cluster(const Cluster &cluster) {
  vector<Cluster> subclusters;
  for (uint c = 0; c < cluster.size(); ++c) {
    const string &name = cluster.get_name(c);
    const string &seq = cluster.get_seq(c);
    size_t i;
    for (i = 0; i < subclusters.size(); i++) {
      float cl = subclusters[i].get_len();
      float sl = seq.size();
      if (min(cl, sl) / max(cl, sl) >= config->min_ratio)
        break;
    }
    if (i == subclusters.size()) {
      subclusters.push_back(
          Cluster(cluster.chrom, cluster.s, cluster.e, cluster.cov));
    }
    subclusters[i].add_seq(name, seq);
  }
  return subclusters;
}

string Caller::run_poa(const vector<string> &seqs) {
  uint n_seqs = seqs.size();
  abpoa_t *ab = abpoa_init();
  abpoa_para_t *abpt = abpoa_init_para();
  abpt->disable_seeding = 1;
  abpt->align_mode = 0; // global
  abpt->out_msa = 0;
  abpt->out_cons = 1;
  abpt->out_gfa = 0;
  // abpt->is_diploid = 1; // TODO: maybe this works now
  abpt->progressive_poa = 0;
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
  string cons = "";
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

/* Call SVs by POA+realignment */
void Caller::pcall(const vector<Cluster> &clusters) {
#pragma omp parallel for num_threads(config->threads) schedule(static, 1)
  for (size_t i = 0; i < clusters.size(); i++) {
    int t = omp_get_thread_num();
    const Cluster &cluster = clusters[i];
    string chrom = cluster.chrom;
    const auto &subclusters = split_cluster(cluster);

    // Sorting clusters by #sequences to get first 2 most weighted clusters
    int i_max1 = -1;
    int i_max2 = -1;
    uint v_max1 = 0;
    uint v_max2 = 0;
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
    vector<int> maxs({i_max1, i_max2});
    // Calling from the two most-weighted clusters
    for (const int i : maxs) {
      if (i == -1)
        continue;
      Cluster cl = subclusters[i];
      if (cl.size() < config->min_cluster_weight)
        continue;

      vector<SV> _svs;

      string ref = string(chromosome_seqs[chrom] + cl.s, cl.e - cl.s + 1);
      string consensus = run_poa(cl.seqs);
      parasail_result_t *result = NULL;
      result = parasail_nw_trace_striped_16(consensus.c_str(), consensus.size(),
                                            ref.c_str(), ref.size(), 10, 1,
                                            &parasail_nuc44);
      parasail_cigar_t *pcigar =
          parasail_result_get_cigar(result, consensus.c_str(), consensus.size(),
                                    ref.c_str(), ref.size(), NULL);
      char *cigar_str = parasail_cigar_decode(pcigar);
      int score = result->score;
      parasail_cigar_free(pcigar);
      parasail_result_free(result);
      _p_alignments[t].push_back(
          Consensus(consensus, cigar_str, chrom, cl.s, cl.e));
      // -- Extracting SVs
      uint rpos = cl.s; // position on reference
      uint cpos = 0;    // position on consensus
      CIGAR cigar;
      cigar.parse_cigar(cigar_str);
      int nv = 0;
      for (const auto cigar_pair : cigar.ops) {
        uint l = cigar_pair.first;
        char op = cigar_pair.second;
        if (op == '=' || op == 'M') {
          rpos += l;
          cpos += l;
        } else if (op == 'I') {
          if (l > config->min_sv_length) {
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
          if (l > config->min_sv_length) {
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
      for (size_t v = 0; v < _svs.size(); v++)
        _svs[v].ngaps = nv;
      for (const SV &sv : _svs)
        _p_svs[t].push_back(sv);
    }
  }
}

/* Merge close and similar SVs */
void Caller::filter_sv_chains() {
  if (svs.size() < 2)
    return;
  sort(svs.begin(), svs.end());

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
  cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
       << endl;
  cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEFAULT"
       << endl;
}
